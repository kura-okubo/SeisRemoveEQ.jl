__precompile__()
module Utils

export get_memoryuse, initlogo, printparams

include("map_removeEQ.jl")
using .Map_removeEQ

using SeisIO, Printf, Dates, JLD2, FileIO



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
        t1 = @elapsed Stest, bt1, bt2 = map_removeEQ(trial_id, InputDict) #[s]
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

    println(mem_per_requestid)
    println(max_num_of_processes_per_parallelcycle)
    println("-------EQ REMOVAL STATS SUMMARY--------")

    println(@sprintf("Number of processes is %d.", nprocs()))

    totaldownloadsize = mem_per_requestid * numofitr
    if totaldownloadsize < MB
        totaldownloadsize = totaldownloadsize * GB / MB #[MB]
        sizeunit = "MB"
    else
        sizeunit = "GB"
    end

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
        println(@sprintf("%24s = %10s", key, string(param["$key"])))
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
