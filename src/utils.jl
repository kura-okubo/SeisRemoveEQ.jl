__precompile__()
module Utils

export convert_tmpfile, defaultinputdict!, printparams, initlogo

using SeisIO, Printf, Dates, JLD2, FileIO


"""
convert_tmpfile(InputDict::Dict)

convert temporal file in "./seisdownload_tmp" to prescribed format.
It has salvage mode, which allows to compile the temporal files in the case of failing during the download.
"""
function convert_tmpfile(InputDict::Dict; salvage::Bool=false)

    println("-------START CONVERTING-------------")

	# save data to fopath file
	fopath = InputDict["fopath"]

	t = jldopen(InputDict["finame"])

	jldopen(fopath, "w") do file
		file["info/DLtimestamplist"] = t["info/DLtimestamplist"];
		file["info/stationlist"]     = t["info/stationlist"];
		file["info/starttime"]       = t["info/starttime"];
		file["info/endtime"]         = t["info/endtime"];
		file["info/DL_time_unit"]    = t["info/DL_time_unit"];
	end

	JLD2.close(t)

	# find all temporal files
    paths = ls(InputDict["tmppath"])
    fmt = InputDict["outputformat"]

    file = jldopen(fopath, "w")

    stationlist     = []
    DLtimestamplist = []
    varnamelist     = []

	IsIsolateComponents = InputDict["IsIsolateComponents"]
	isostationlist = []

    for path in paths
        #println(path)

		# split path
		tmpname = split(path, "/")[end]
		y, d, net, sta, loc, cha, _, _ = split(tmpname, ".")
		iso_stationinfo = (join([y, d, net, sta, loc], "-"), cha)
		#println(iso_stationinfo)

        try
            S = SeisIO.rseis(path)[1]
            #println(S)

            for ii = 1:S.n #loop at each seis channel

				# check channel isolation
				conflictsta = filter(x -> x[1]==iso_stationinfo[1] && string(x[2][end][end])==string(iso_stationinfo[2][end]), isostationlist)
				#println(conflictsta)
				if !isempty(conflictsta)
					# here this channel has conflicting channel such as same day, same components but different channel.
					# if IsIsolateComponents == true, skip to avoid multiple stations for the purpose of cross-correlation.
					if IsIsolateComponents
						sta1 = join([iso_stationinfo[1], iso_stationinfo[2]], "-")
						sta2 = join([conflictsta[1][1], conflictsta[1][2]], "-")
						txt = @sprintf("Two identical location but different channel were found.\n%s :%s: Discard %s.",
									sta1, sta2, sta1)
						println(txt)
						#rm(path)
						continue;
					end
				end

				# update isostationlist
				push!(isostationlist, iso_stationinfo)

                # make station list
                staid = S[ii].id

                # save data (checking whether it's already in the jld2 because it causes an error)
                #parse info
                s_str = string(u2d(S[ii].t[1,2]*1e-6))

                # select output format
                if fmt == "JLD2"
                    yj = parse(Int64, s_str[1:4])
                    mj = parse(Int64, s_str[6:7])
                    dj = parse(Int64, s_str[9:10])
                    tj = string(s_str)[11:19]

                    djm2j = md2j(yj, mj, dj)
                    groupname = string(yj)*"."*string(djm2j)*"."*tj #Year_Julianday_Starttime
                    varname = joinpath(groupname, staid)

                    if isempty(filter(x -> x==varname, varnamelist))
                        push!(varnamelist, varname)
                        file[varname] = S[ii]
                    end
                    # @info "save data $varname"
                else
                    error("output format in $fmt is not implemented yet.")
                end
            end

			rm(path)

        catch y
            println(y)
        end
    end

    JLD2.close(file)

	rm(InputDict["tmppath"], recursive=true, force=true)

    return nothing
end

"""
printparams(param::Dict)

print parameters
"""
function printparams(param::Dict)
    printstyled("---Input Parameters---\n"; color=:cyan, bold=true)
    for key in keys(param)
        println(@sprintf("%-24s = %-10s", key, string(param["$key"])))
    end
end

"""
defaultinputdict(InputDict::Dict)

default if input parameter is missing.
"""
function defaultinputdict!(InputDict::Dict)

	# default values
	def = Dict()
	def["finame"]				 	= "./input.jld2"
	def["IsKurtosisRemoval"] 		= true
	def["max_edgetaper_duration"] 	= 60 * 5
	def["kurtosis_tw_sparse"] 		= 60
	def["kurtosis_timewindow"] 		= 60*3
	def["kurtosis_threshold"] 		= 2.0
	def["IsSTALTARemoval"] 			= true
	def["stalta_longtimewindow"] 	= 60*60*2
	def["stalta_threshold"] 		= 1.2
	def["stalta_absoluteclip"] 		= 1e20
	def["max_wintaper_duration"] 	= 60 * 3
	def["removal_shorttimewindow"] 	= 60 * 3
	def["overlap"] 					= 60
	def["IsIsolateComponents"] 		= false
	def["IsSaveFig"] 				= false
	def["plot_kurtosis_α"] 			= 1.2
	def["plot_boxheight	"] 			= 1.5
	def["plot_span"] 				= 100
	def["outputformat"]				= "JLD2"
	def["IsStartendtime"] 			= false
	def["fodir"] 					= "./dataset"
    def["foname"] 					= "eq_removed.jld2"

	for key in keys(def)
		if !haskey(InputDict, key)
			InputDict["$key"] = def["$key"]
		end
	end

end



"""
initlogo()

print initial logo
"""
function initlogo()

    print("
    _____      _      ______
   /  ___|    (_)     | ___ \\
   \\ `--.  ___ _ ___  | |_/ /___ _ __ ___   _____   _____
    `--. \\/ _ \\ / __| |    // _ \\ '_ ` _ \\ / _ \\ \\ / / _ \\
   /\\__/ /  __/ \\__ \\ | |\\ \\  __/ | | | | | (_) \\ V /  __/
   \\____/ \\___|_|___/ \\_| \\_\\___|_| |_| |_|\\___/ \\_/ \\___|
    ______ ____
   |  ____/ __ \\       |
   | |__ | |  | |      |
   |  __|| |  | |      |  v1.0 (Last update 07/07/2019)
   | |___| |__| |      |  © Kurama Okubo
   |______\\___\\_\\      |

")

    println("Job start running at "*string(now())*"\n")

end

end
