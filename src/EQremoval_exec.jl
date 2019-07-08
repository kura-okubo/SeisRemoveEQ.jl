#using Distributed
#addprocs(3)

using Dates, JLD2

@everywhere include("./lib/ParallelEQremoval.jl")
@everywhere using SeisIO, .EQremoval

#==================================================#

finame      = "./dataset/BPnetwork.jld2" # data is saved at ./dataset/$foname.jld2

max_num_of_processes_per_parallelcycle = 10 #define how many time stamps you parallelization within one pmap cycle: reduce this to reduce memory allocation

maxedgetaperduration = 60 * 5 #[s] tapering to avoid edge effect

kurtosis_timewindow = 60 * 5 #300 # kurtosis timewindow[s]
threshold = 2.0 #kurtosis threshold
overlap = 30; # overlap of kurtosis removal window[s]

#invert_tukey window
tukey_α = 0.1

#---Plotting variables---#
kurtosis_α = 0.5 # amplitude normalization for plotting kurtosis
boxheight = 1.5

span = 100 # to reduce number of plot point

foutputname = "./dataset/BPnetwork_removedEQ.jld2"

IsSaveFig = true # save figure to check removal
#----------------------------------------------------#

#==================================================#

t = jldopen(finame)

#save data to foutput jld2 file
jldopen(foutputname, "w") do file
        file["info/DLtimestamplist"] = t["info/DLtimestamplist"];
        file["info/stationlist"]     = t["info/stationlist"];
        file["info/starttime"]       = t["info/starttime"];
        file["info/endtime"]         = t["info/endtime"];
        file["info/DL_time_unit"]    = t["info/DL_time_unit"];
end

#intended to parallelization
InputDictionary = Dict( "finame"                => finame,
                        "maxedgetaperduration"  => maxedgetaperduration,
                        "kurtosis_timewindow"   => kurtosis_timewindow,
                        "threshold"             => threshold,
                        "overlap"               => overlap,
                        "tukey_α"               => tukey_α,
                        "kurtosis_α"            => kurtosis_α,
                        "boxheight"             => boxheight,
                        "span"                  => span,
                        "IsSaveFig"             => IsSaveFig)



DLtimestamplist = t["info/DLtimestamplist"];
stationlist     = t["info/stationlist"];
NumofTimestamp  = length(DLtimestamplist)

JLD2.close(t)



file = jldopen(foutputname, "r+")

E = []

bt_getkurtosis = 0.0
bt_removeeq = 0.0

if max_num_of_processes_per_parallelcycle >= NumofTimestamp

		EE = pmap(x -> ParallelEQremoval(x, InputDictionary), 1:NumofTimestamp)

		for i = 1:size(EE)[1]
			push!(E, EE[i][1])
			global bt_getkurtosis += EE[i][2]
		    global bt_removeeq += EE[i][3]
		end

		# save data to jld2
	   for ii = 1:size(E)[1] #loop at each starttime
		   for jj = 1:size(E[1])[1] #loop at each station id

			   requeststr =E[ii][jj].id
			   varname = joinpath(DLtimestamplist[ii], requeststr)
			   #save_SeisData2JLD2(fopath, varname, S[ii][jj])
			   file[varname] = E[ii][jj]

		   end
	   end

else

	#parallelization by time
	pitr = 1

	while pitr <=  NumofTimestamp

	    startid1 = pitr
	    startid2 = pitr + max_num_of_processes_per_parallelcycle - 1

	    if startid1 == NumofTimestamp
	        #use one
	        E, bt1, bt2 = pmap(x -> ParallelEQremoval(x, InputDictionary), startid1:startid1)

	    elseif startid2 <= NumofTimestamp
	        #use all processors
	        E, bt1, bt2 = pmap(x -> ParallelEQremoval(x, InputDictionary), startid1:startid2)

	    else
	        #use part of processors
	        startid2 = startid1 + mod(NumofTimestamp, nprocs()) - 1
	        println(startid2)
	        E, bt1, bt2 = pmap(x -> ParallelEQremoval(x, InputDictionary), startid1:startid2)
	    end
	    # save data to jld2
	    for ii = 1:size(E)[1] #loop at each starttime
	        for jj = 1:size(E[1])[1] #loop at each station id

	            requeststr =E[ii][jj].id
	            varname = joinpath(DLtimestamplist[startid1+ii-1], requeststr)
	            #save_SeisData2JLD2(fopath, varname, S[ii][jj])
	            file[varname] = E[ii][jj]

	        end
	    end

	    global pitr += max_num_of_processes_per_parallelcycle
	end
end

JLD2.close(file)

println("bt_getkurtosis=$(bt_getkurtosis)[s]")
println("bt_removeEQ=$(bt_removeeq)[s]")

print("process done at:")
println(now())
str = "EQ is successfully removed from raw data: output = $foutputname"
printstyled(str; color=:green, bold=true)
