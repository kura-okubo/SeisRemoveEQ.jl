__precompile__()
module Utils

export get_memoryuse, initlogo, printparams

include("map_removeEQ.jl")
using .Map_removeEQ

using SeisIO, Printf, Dates, JLD2, FileIO, Distributed



const KB = 1024.0 #[bytes]
const MB = 1024.0 * KB
const GB = 1024.0 * MB

"""
get_memoryuse(InputDict::Dict)

evaluate memory use and estimate computational time
"""
function get_memoryuse(InputDict::Dict)

    trial_id        = 1

    while true
        global t1 = @elapsed EE = map_removeEQ(trial_id, InputDict) #[s]

		global Stest = []

        for i = 1:size(EE[1])[1]
			push!(Stest, EE[1][i])
		end

        dl = [Stest[i].misc["dlerror"] for i in 1:size(Stest)[1]]
        if issubset(0, dl)
            break;
        else
            trial_id += 1
        end
    end

    numofitr = InputDict["NumofTimestamp"] # number of iteration by parallel processing
    if trial_id == numofitr - 1
        error("all requests you submitted with input dictionary was failed. Please check the station availability in your request.")
    end

    mem_per_requestid = 1.2 * sizeof(Stest) / GB #[GB] *for the safty, required memory is multiplied by 1.2

    max_num_of_processes_per_parallelcycle = floor(Int64, InputDict["MAX_MEM_PER_CPU"]/mem_per_requestid)
    estimated_downloadtime = now() + Second(round(3 * t1 * numofitr / nprocs()))

	printstyled("---EQ REMOVAL STATS SUMMARY---\n"; color=:cyan, bold=true)

    println(@sprintf("Number of processes is %d.", nprocs()))

    totaldownloadsize = mem_per_requestid * numofitr
    if totaldownloadsize < MB
        totaldownloadsize = totaldownloadsize * GB / MB #[MB]
        sizeunit = "MB"
    else
        sizeunit = "GB"
    end

	println(@sprintf("One process returns %4.2e [%s] of dataset so maximum num of parallel processes is %d.",
	 			totaldownloadsize, sizeunit, max_num_of_processes_per_parallelcycle))
	println(@sprintf("Download will finish at %s.", round(estimated_downloadtime, Dates.Second(1))))
    println("*We have a time lag with processing time above, like in 10 minutes or so.*")
    println("*This estimation also changes if some processes fail and are skipped.*")

    return max_num_of_processes_per_parallelcycle

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
	def["IsSaveFig"] 				= false
	def["plot_kurtosis_α"] 			= 1.2
	def["plot_boxheight	"] 			= 1.5
	def["plot_span"] 				= 100
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
