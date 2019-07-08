module Map_removeEQ

export map_removeEQ

using Dates, JLD2, DSP, PlotlyJS, Printf, ORCA, FileIO

using SeisIO

include("get_kurtosis.jl")
include("remove_eq.jl")
using .Get_kurtosis, .Remove_eq

"""
    ParallelEQremoval(dlid, InputDict::Dict)
    remove earthquake and save it into jld2 file.

"""
function map_removeEQ(dlid, InputDict::Dict)

    #store data
    finame                 =InputDict["finame"]
    max_edgetaper_duration =InputDict["max_edgetaper_duration"]
    kurtosis_timewindow    =InputDict["kurtosis_timewindow"]
    kurtosis_threshold     =InputDict["kurtosis_threshold"]
    overlap                =InputDict["overlap"]
    invert_tukey_α         =InputDict["invert_tukey_α"]
    plot_kurtosis_α        =InputDict["plot_kurtosis_α"]
    plot_boxheight         =InputDict["plot_boxheight"]
    plot_span              =InputDict["plot_span"]
    fodir                  =InputDict["fodir"]
    foname                 =InputDict["foname"]
    fopath                 =InputDict["fopath"]
    IsSaveFig              =InputDict["IsSaveFig"]

    DLtimestamplist        =InputDict["DLtimestamplist"]
    stationlist            =InputDict["stationlist"]
    NumofTimestamp         =InputDict["NumofTimestamp"]

    tstamp = DLtimestamplist[dlid]

    if mod(dlid, round(0.1*NumofTimestamp)+1) == 0
        println(@sprintf("start process %s", tstamp))
    end

    SRall = SeisData(length(stationlist))

    icount = 0

    bt_getkurtosis = 0.0
    bt_removeeq = 0.0

    for st = stationlist
        #S = t[joinpath(tstamp, st)]
        S = FileIO.load(finame, joinpath(tstamp, st))

        if S.misc["dlerror"] == 0
            dt = 1/S.fs
            tvec = collect(0:S.t[2,1]-1) * dt ./ 60 ./ 60

            #tapering to avoid instrumental edge artifacts
            SeisIO.taper!(S,  t_max = max_edgetaper_duration, α=0.05)
            S1 = deepcopy(S)

            bt_1 = @elapsed S1 = Get_kurtosis.get_kurtosis(S1, float(kurtosis_timewindow))
            tw = kurtosis_timewindow; #so far this is most stable
            bt_2 = @elapsed S1 = Remove_eq.detect_eq(S1, float(tw), float(kurtosis_threshold), float(overlap))
            bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, S, float(invert_tukey_α), plot_kurtosis_α,
                                plot_boxheight, plot_span, fodir, tstamp, tvec, IsSaveFig)

            #remove kurtosis for memory
            S1.misc["kurtosis"] = []

            bt_getkurtosis += bt_1
            bt_removeeq += bt_2 + bt_3

        else
            #download error found: save as it is.
            S1 = S
        end

        icount += 1
        SRall[icount] = S1

    end

    return (SRall, bt_getkurtosis, bt_removeeq)
end

end
