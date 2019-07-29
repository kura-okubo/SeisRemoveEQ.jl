module Map_removeEQ

export map_removeEQ

using Dates, JLD2, DSP, PlotlyJS, Printf, ORCA, FileIO, Suppressor

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
    IsKurtosisRemoval      =InputDict["IsKurtosisRemoval"]
    max_edgetaper_duration =InputDict["max_edgetaper_duration"]
    kurtosis_tw_sparse     =InputDict["kurtosis_tw_sparse"]
    kurtosis_timewindow    =InputDict["kurtosis_timewindow"]
    kurtosis_threshold     =InputDict["kurtosis_threshold"]
    IsSTALTARemoval        =InputDict["IsSTALTARemoval"]
    stalta_longtimewindow  =InputDict["stalta_longtimewindow"]
    stalta_threshold       =InputDict["stalta_threshold"]
    stalta_absoluteclip    =InputDict["stalta_absoluteclip"]
    max_wintaper_duration  =InputDict["max_wintaper_duration"]
    removal_shorttimewindow=InputDict["removal_shorttimewindow"]
    overlap                =InputDict["overlap"]
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

    SRall = []

    icount = 0

    bt_getkurtosis = 0.0
    bt_removeeq = 0.0

    for st = stationlist
        #S = t[joinpath(tstamp, st)]
        st1 = replace(st, "-"=>"")

        loaddata(finame, path) = try
            @suppress return FileIO.load(finame, path)
        catch
           return false;
        end

        # load raw data
        Stemp = loaddata(finame, joinpath(tstamp, st1))

        # skip if it does not exist
        if Stemp == false
            continue
        end

        #convert if it's SeisChannel

        if Stemp isa SeisChannel
            Stemp = SeisData(Stemp)
        end

        SRtemp = SeisData(size(Stemp)[1])

        #loop by SeisChannel
        for Seisid = 1:size(Stemp)[1]

            S = Stemp[Seisid]

            if !haskey(S.misc, "dlerror")
                if iszero(S.x)
                    S.misc["dlerror"] = true
                else
                    S.misc["dlerror"] = false
                end
            end

            if S.misc["dlerror"] == 0
                dt = 1/S.fs
                #tvec = collect(0:S.t[2,1]-1) * dt ./ 60 ./ 60
                tvec = collect(0:length(S.x)-1) * dt ./ 60 ./ 60

                #tapering to avoid instrumental edge artifacts
                SeisIO.taper!(S,  t_max = max_edgetaper_duration, α=0.05)
                S1 = deepcopy(S)

                #set long window length to user input since last window of previous channel will have been adjusted
                S1.misc["eqtimewindow"] = fill(true, length(S1.x))

                if IsKurtosisRemoval
                    # compute kurtosis and detect earthqukes
                    bt_1 = @elapsed S1 = Get_kurtosis.get_kurtosis(S1, float(kurtosis_timewindow), float(kurtosis_tw_sparse))
                    bt_2 = @elapsed S1 = Remove_eq.detect_eq_kurtosis(S1, tw=float(removal_shorttimewindow), kurtosis_threshold=float(kurtosis_threshold), overlap=float(overlap))

                    btsta_1 = 0

                    if IsSTALTARemoval
                        # detect earthquake and tremors by STA/LTA
                        btsta_1 = @elapsed S1 = Remove_eq.detect_eq_stalta(S1, float(stalta_longtimewindow), float(removal_shorttimewindow),
                                            float(stalta_threshold), float(overlap), float(stalta_absoluteclip))
                    end

                    bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, S, float(plot_kurtosis_α), float(max_wintaper_duration),
                                    plot_boxheight, trunc(Int, plot_span), fodir, tstamp, tvec, IsSaveFig)


                    bt_getkurtosis += bt_1
                    bt_removeeq += bt_2 + bt_3 + btsta_1

                    #println([bt_2, bt_3, btsta_1])
                    #if mod(dlid, round(0.1*NumofTimestamp)+1) == 0
                    #    println([bt_1, bt_2, btsta_1, bt_3])
                    #end

                    #remove kurtosis for reduce size
                    S1.misc["kurtosis"] = []
                    S1.misc["eqtimewindow"] = []

                else
                    #only STA/LTA
                    if IsSTALTARemoval

                        bt_2 = @elapsed S1 = Remove_eq.detect_eq_stalta(S1, float(stalta_longtimewindow), float(removal_shorttimewindow),
                                            float(stalta_threshold), float(overlap))
                        bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, S, float(plot_kurtosis_α), float(max_wintaper_duration),
                                        plot_boxheight, trunc(Int, plot_span), fodir, tstamp, tvec, IsSaveFig)

                        #remove kurtosis for reduce size
                        S1.misc["kurtosis"] = []
                        S1.misc["eqtimewindow"] = []

                        bt_getkurtosis += 0.0
                        bt_removeeq += bt_2 + bt_3

                    else
                        #no removal process.
                        println("Both 'IsKurtosisRemoval' and 'IsKurtosisRemoval' are false. No removal process is executed.")
                        exit(0)
                    end

                end

            else
                #download error found: save as it is.
                S1 = S
            end

            SRtemp[Seisid] = S1

        end

        icount += 1

        push!(SRall, SRtemp)

    end

    return (SRall, bt_getkurtosis, bt_removeeq)
end

end
