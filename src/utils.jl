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

	file = jldopen(fopath, "w")

	#!!!should be debuged because of Isocomponents!!!#
	file["info/stationlist"]     = t["info/stationlist"];

	if InputDict["IsStartendtime"]
		file["info/DLtimestamplist"] = InputDict["DLtimestamplist_selected"];
		file["info/starttime"]       = InputDict["starttime"];
		file["info/endtime"]         = InputDict["endtime"];
	else
		file["info/DLtimestamplist"] = t["info/DLtimestamplist"];
		file["info/starttime"]       = t["info/starttime"];
		file["info/endtime"]         = t["info/endtime"];
	end

	JLD2.close(t)

	# find all temporal files
    paths_all = ls(InputDict["tmppath"])
    fmt = InputDict["outputformat"]

    stationlist     = []
    DLtimestamplist = []
    varnamelist     = []

	if InputDict["IsIsolateComponents"]
		# isolate components based on priority dictionary
		paths = isolate_components(paths_all, InputDict)
	else
		paths = paths_all
	end

    for path in paths
        #println(path)

		# split path
		# tmpname = split(path, "/")[end]

		# ftmpname = split(tmpname, ".")
		# if occursin("-", ftmpname[3])
		# 	# format would be y, jd, T00:00:00, sta, loc, cha
		# 	y, d, tmpT, net, sta, loc, cha = split(tmpname, ".")
		# 	#iso_stationinfo = (join([y, d, net, sta, loc], "-"), cha)
		# 	iso_stationinfo = (join([y, d, net, sta], "-"), cha)
		#
		# elseif !occursin("-", ftmpname[3]) && IsIsolateComponents == true
		# 	@warn "format of tmp file is not y, jd, time, sta, loc, cha. Turn off IsIsolateComponents."
		# 	IsIsolateComponents = false
		# 	iso_stationinfo = []
		# end

		#println(iso_stationinfo)


        S = try
				SeisIO.rseis(path)[1]
			catch y
				println(y)
				@warn("cannot read tmpfile in seisremoveeq_tmp_sample. skipping")
				continue;
			end
        #println(S)

        for ii = 1:S.n #loop at each seis channel

			# if IsIsolateComponents
			# 	# check channel isolation: if y, d, net, sta and components are same, we discard one of them.
			# 	@show isostationlist
			# 	@show iso_stationinfo
			# 	conflictsta = findfirst(x -> x[1]==iso_stationinfo[1] && string(x[2][end])==string(iso_stationinfo[2][end]), isostationlist)
			# else
			# 	conflictsta = []
			# end

			#println(conflictsta)
			# if !isempty(conflictsta)
			# 	# here this channel has conflicting channel such as same day, same components but different channel.
			# 	# if IsIsolateComponents == true, skip to avoid multiple stations for the purpose of cross-correlation.
			# 	if IsIsolateComponents
			# 		sta1 = join([iso_stationinfo[1], iso_stationinfo[2]], "-")
			# 		sta2 = join([conflictsta[1][1], conflictsta[1][2]], "-")
			# 		txt = @sprintf("Two identical location but different channel were found.\n%s :%s: Discard %s.",
			# 					sta1, sta2, sta1)
			# 		println(txt)
			# 		#rm(path)
			# 		continue;
			# 	end
			# end

			# if IsIsolateComponents
			# 	# update isostationlist
			# 	push!(isostationlist, iso_stationinfo)
			# end

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

            else
                error("output format in $fmt is not implemented yet.")
            end
        end

		# remove tmpfile
		rm(path)

    end

    JLD2.close(file)

	rm(InputDict["tmppath"], recursive=true, force=true)

    return nothing
end

"""
isolate_components(paths_all::AbstractArray, InputDict::Dict)

Isolate components based on InputDict["priority_channles"].

I. What is the purpose?
	Some stations have multiple stations with different channels
	(e.g. BP.LCCB..BP1 and BP.LCCB.40.SP1 during transition.)
	This increases number of xcorr processes, while the phisicall
	meaning is auto correlation using different stations; thus we duplicate
	auto-correlation in this case (e.g. BP1-BP1 and SP1-BP1 at same place).
	To avoid that, we implemented this function to isolate components with
	each stations based on priority, so that we can perform long-term dv/v
	analysis without such duplications.

II. Process flow

	1. Try to find channels which has same network.station, and same component
	but different channel

	2. If found, search priority with `haskey(InputDict["priority_channles"])`.
	If not found, arbitrually pick up the first one.

	3. If the second one has priority, search the order of it in Dictionary;
	(e.g. findfirst("BP1", InputDict["priority_channles"]["BP"])). if not,
	arbitrually pick up the first one.

	4. Then search the priority for first one;
	(e.g. findfirst("SP1", InputDict["priority_channles"]["BP"])).
	If not found, arbitrually pick up the second one as it is listed in
	priority dictionary.

	5. Compare the priority between first and second one, and pick up the
	earlier one as representive channel at this station.

III. Potential issue

	We assume that the physical measurement at same station does not change so much
	that the cross-correlation is not influenced by the replacement of channels.
	If it is not satisfied, it causes diference in the cross-correlatino result.
	Please check the consistency between channels at same station if you find
	some discontinuous result before and after switch of channel.

"""
function isolate_components(paths_all::AbstractArray, InputDict::Dict)

	iso_list = []

	for path in paths_all
		tmp = split(path, "/")[end]
		# read meta data from file name
		ftmpname = split(tmp, ".")

		if occursin("-", ftmpname[3])
			# format would be y, jd, T00-00-00, sta, loc, cha
			y, d, tmpT, net, sta, loc, cha = split(ftmpname, ".")
			#iso_stationinfo = (join([y, d, net, sta, loc], "-"), cha)
			# (time, net, sta for find the station, channel name and component)
			push!(iso_list, (path, join([y, d, tmpT, net, sta], "-"), net, sta, cha[1:2], cha[3]))

		else
			@warn "Format of tmp file is not y, jd, time, sta, loc, cha. Turn off IsIsolateComponents."
			return paths
		end
	end

	@show iso_list

	isocomp_idlist = []

	for (ista, current_sta) in enumerate(iso_list)
		for jsta = ista:length(iso_list)
			compared_sta = iso_list[j]
			if current_sta[2] != compared_sta[2] || current_sta[6] != compared_sta[6]
				# there is no conflict in channel, so add the current one to isocomp list
				if ista ∉ isocomp_idlist
					push!(isocomp_idlist, ista)
				end
			else
				# here current_sta[2] == compared_sta[2] && current_sta[6] == compared_sta[6]
				# perform process 2. we currently take into account priority with
				# all stations in certain network

				iso_net = compared_sta[3]
				current_cha = current_sta[5]
				compared_cha = compared_sta[5]

				if !haskey(InputDict["priority_channles"],iso_net)
					# second compared station does not have priority. take the current one
					if ista ∉ isocomp_idlist
						push!(isocomp_idlist, ista)
					end

				else
					# this has priority; perform process 3. e.g. compared_sta[5] = "SP"
					priority_compared = findfirst(compared_cha, InputDict["priority_channles"][iso_net])
					if isempty(priority_compared)
						# this channel has no priority. take the current one
						if ista ∉ isocomp_idlist
							push!(isocomp_idlist, ista)
						end

					else
						# this has priority so that search priority for current one.
						priority_current = findfirst(current_cha, InputDict["priority_channles"][iso_net])
						if isempty(priority_current)
							# current channel has no priority. take the compared one
							if jsta ∉ isocomp_idlist
								push!(isocomp_idlist, jsta)
							end

							# filter out ista if it's in iscomp_list
							filter!(x -> x != ista, isocomp_idlist)

						else
							# both current and compared one has priority. compare the order, and
							# add or replace ista in isocomp_idlist
							if priority_current > priority_compared
								# take current one
								if ista ∉ isocomp_idlist
									push!(isocomp_idlist, ista)
								end
							else
								# take compared one
								if jsta ∉ isocomp_idlist
									push!(isocomp_idlist, jsta)
								end
								# filter out ista if it's in iscomp_list
								filter!(x -> x != ista, isocomp_idlist)
							end
						end
					end
				end
			end
		end
	end

	#DEBUG:
	for id in isocomp_idlist
		temppath = paths_all[id]
		tmp = split(temppath, "/")[end]
		println(tmp)
	end

	return paths_all[isocomp_idlist]

end

"""
printparams(param::Dict)

print parameters
"""
function printparams(param::Dict)
    printstyled("-----------Input Parameters-----------\n"; color=:cyan, bold=true)
    for key in keys(param)

		try
	        if length(string(param["$key"])) > 60
	            param_str = string(param["$key"])[1:30]*"..."
	        else
	            param_str = string(param["$key"])
	        end
	        println(@sprintf("%-24s = %-10s", key, param_str))
		catch
		end
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
	def["priority_channles"] 		= Dict()
	def["IsSaveFig"] 				= false
	def["plot_kurtosis_α"] 			= 1.2
	def["plot_boxheight	"] 			= 1.5
	def["plot_span"] 				= 100
	def["outputformat"]				= "JLD2"
	def["IsStartendtime"] 			= false
	def["fodir"] 					= "./dataset"
	def["foname"] 					= "eq_removed.jld2"

	#For dump everytrace
	#This is always false. If you need to plot some figures for
	#details of kurtosis and STA/LTA, please turn it true; it slowdowns down computation.
	def["dumptraces"] 				= false
	def["dumppath"]					= InputDict["fodir"]*"/dumptraces"

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
