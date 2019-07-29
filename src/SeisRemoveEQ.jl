__precompile__()
module SeisRemoveEQ

include("utils.jl")
include("map_removeEQ.jl")

using .Utils
using .Map_removeEQ

using Dates, Printf, JLD2, FileIO, Distributed

export seisremoveEQ

function seisremoveEQ(InputDict::Dict)

	Utils.initlogo()

	Utils.defaultinputdict!(InputDict)

	finame = InputDict["finame"]
	fodir  = InputDict["fodir"]
	foname = InputDict["foname"]

	mkpath(fodir)
	fopath=joinpath(fodir, foname*".jld2")
	InputDict["fopath"] = fopath

	#save data to fopath file
	t = jldopen(finame)
	jldopen(fopath, "w") do file
		file["info/DLtimestamplist"] = t["info/DLtimestamplist"];
		file["info/stationlist"]     = t["info/stationlist"];
		file["info/starttime"]       = t["info/starttime"];
		file["info/endtime"]         = t["info/endtime"];
		file["info/DL_time_unit"]    = t["info/DL_time_unit"];
	end

	DLtimestamplist = t["info/DLtimestamplist"];
	stationlist     = t["info/stationlist"];
	NumofTimestamp  = length(DLtimestamplist)
	JLD2.close(t)

	#print parameters
	printparams(InputDict)

	InputDict["DLtimestamplist"] = DLtimestamplist
	InputDict["stationlist"] = stationlist
	InputDict["NumofTimestamp"] = NumofTimestamp

	#evaluate memory use
	printstyled("---Start Test Process---\n"; color=:cyan, bold=true)

	max_num_of_processes_per_parallelcycle = get_memoryuse(InputDict)

	if max_num_of_processes_per_parallelcycle < 1
		error("Memory allocation is not enought (currently $MAX_MEM_PER_CPU [GB]). Please inclease MAX_MEM_PER_CPU or decrease number of stations")
	end

	printstyled("---Start removing EQ---\n"; color=:cyan, bold=true)

	# parallelize process by time stamp

	file = jldopen(fopath, "r+")

	E = []
	bt_getkurtosis = 0.0
	bt_removeeq = 0.0

	if max_num_of_processes_per_parallelcycle >= NumofTimestamp
		# one process cycle cover all time stamp
		EE = pmap(x -> map_removeEQ(x, InputDict), 1:NumofTimestamp)

		for i = 1:size(EE)[1]
			push!(E, EE[i][1])
			bt_getkurtosis += EE[i][2]
		    bt_removeeq += EE[i][3]
		end

		# save data to jld2
		for ii = 1:length(E) #loop at each starttime
			for jj = 1:length(E[ii]) #loop at each station id
			   requeststr =E[ii][jj][1].id
			   varname = joinpath(DLtimestamplist[ii], requeststr)
			   #save_SeisData2JLD2(fopath, varname, S[ii][jj])
			   file[varname] = E[ii][jj] #SeisData
			end
		end

	else
		# processes by each max_num_of_processes_per_parallelcycle
		pitr = 1

		while pitr <=  NumofTimestamp

			E = []

		    startid1 = pitr
		    startid2 = pitr + max_num_of_processes_per_parallelcycle - 1

		    if startid1 == NumofTimestamp
		        #use one
		        EE = pmap(x -> map_removeEQ(x, InputDict), startid1:startid1)

		    elseif startid2 <= NumofTimestamp
		        #use all processors
		        EE = pmap(x -> map_removeEQ(x, InputDict), startid1:startid2)

		    else
		        #use part of processors
		        startid2 = startid1 + mod(NumofTimestamp, max_num_of_processes_per_parallelcycle) - 1
		        println(startid2)
		        EE = pmap(x -> map_removeEQ(x, InputDict), startid1:startid2)
		    end

			println(EE[1][1])

			for i = 1:size(EE)[1]
				push!(E, EE[i][1])
				bt_getkurtosis += EE[i][2]
			    bt_removeeq += EE[i][3]
			end

		    # save data to jld2

			# save data to jld2
			for ii = 1:length(E) #loop at each starttime
				for jj = 1:length(E[ii]) #loop at each station id
				   requeststr =E[ii][jj][1].id
				   varname = joinpath(DLtimestamplist[ii], requeststr)
				   #save_SeisData2JLD2(fopath, varname, S[ii][jj])
				   file[varname] = E[ii][jj] #SeisData
				end
			end

		    pitr += max_num_of_processes_per_parallelcycle

		end
	end

	JLD2.close(file)

	printstyled("---Summary---\n"; color=:cyan, bold=true)
	println("time to get kurtosis =$(bt_getkurtosis)[s]")
	println("time to remove EQ    =$(bt_removeeq)[s]")
	print("\nprocess done at:")
	println(now())
	str = "EQ is successfully removed from raw data:\nOutput = $fopath\n"
	printstyled(str; color=:green, bold=true)

end

end # module
