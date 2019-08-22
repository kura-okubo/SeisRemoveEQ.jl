module Map_removeEQ

export map_removeEQ

using SeisIO, Dates, JLD2, PlotlyJS, Printf, FileIO

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
    plot_fmt               =InputDict["plot_fmt"]
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

    bt_getkurtosis = 0.0
    bt_removeeq = 0.0

    for st = stationlist
        #S = t[joinpath(tstamp, st)]
        st1 = replace(st, "-"=>"")

        loaddata(finame, path) = try
            S = SeisData()
            S = FileIO.load(finame, path)
            return S
        catch
            return false;
        end

        # load raw data
        Stemp = loaddata(finame, joinpath(tstamp, st1))

        # skip if it does not exist
        if Stemp == false
            continue;
        end

        #convert if it's SeisChannel
        if Stemp isa SeisIO.SeisChannel
            Stemp = SeisData(Stemp)
        end

        SremEQ = SeisData()

        #loop by SeisChannel and save into temporal wile
        for Seisid = 1:Stemp.n

            Schan = Stemp[Seisid]

            if !haskey(Schan.misc, "dlerror")
                if iszero(Schan.x)
                    Schan.misc["dlerror"] = true
                else
                    Schan.misc["dlerror"] = false
                end
            end

            if Schan.misc["dlerror"] == 0

                dt = 1/Schan.fs
                tvec = collect(0:length(Schan.x)-1) * dt ./ 60 ./ 60

                #tapering to avoid instrumental edge artifacts
                SeisIO.taper!(Schan,  t_max = max_edgetaper_duration, α=0.05)
                S1 = deepcopy(Schan)

                #set long window length to user input since last window of previous channel will have been adjusted
                S1.misc["eqtimewindow"] = fill(true, length(S1.x))

                if IsKurtosisRemoval
                    # compute kurtosis and detect
                    bt_1 = @elapsed S1 = Get_kurtosis.get_kurtosis(S1, float(kurtosis_timewindow), float(kurtosis_tw_sparse))

                    bt_2 = @elapsed S1 = Remove_eq.detect_eq_kurtosis(S1, tw=float(removal_shorttimewindow), kurtosis_threshold=float(kurtosis_threshold), overlap=float(overlap))

                    btsta_1 = 0

                    if IsSTALTARemoval
                        # detect earthquake and tremors by STA/LTA
                        btsta_1 = @elapsed S1 = Remove_eq.detect_eq_stalta(S1, float(stalta_longtimewindow), float(removal_shorttimewindow),
                                            float(stalta_threshold), float(overlap), float(stalta_absoluteclip))
                    end

                    bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, Schan, float(plot_kurtosis_α), float(max_wintaper_duration),
                                    plot_boxheight, trunc(Int, plot_span), plot_fmt, fodir, tstamp, tvec, IsSaveFig)

                    bt_getkurtosis += bt_1
                    bt_removeeq += bt_2 + bt_3 + btsta_1

                    #println([bt_2, bt_3, btsta_1])
                    #if mod(dlid, round(0.1*NumofTimestamp)+1) == 0
                    #    println([bt_1, bt_2, btsta_1, bt_3])
                    #end

                else
                    #only STA/LTA
                    if IsSTALTARemoval

                        data.misc["kurtosis"] = zeros(length(S1.x))

                        bt_2 = @elapsed S1 = Remove_eq.detect_eq_stalta(S1, float(stalta_longtimewindow), float(removal_shorttimewindow),
                                            float(stalta_threshold), float(overlap))
                        bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, S, float(plot_kurtosis_α), float(max_wintaper_duration),
                                        plot_boxheight, trunc(Int, plot_span), plot_fmt, fodir, tstamp, tvec, IsSaveFig)


                        bt_getkurtosis += 0.0
                        bt_removeeq += bt_2 + bt_3

                    else
                        #no removal process.
                        @warn "Both 'IsKurtosisRemoval' and 'IsKurtosisRemoval' are false. No removal process is executed. Abort."
                        exit(0)
                    end

                end

                #it's not allowed to save this into binary;
                delete!(S1.misc, "eqtimewindow")
                delete!(S1.misc, "kurt")

            else
                #download error found: save as it is.
                S1 = Schan
            end

            SremEQ += S1

        end


        # if some of SeisChannels in Stemp have a data, save temp file
        y, jd = parse.(Int64, split(InputDict["DLtimestamplist"][dlid], "."))
        fname_out = join([tstamp, st1, "FDSNWS","dat"], '.')
        # save as intermediate binary file
        t_write = @elapsed SeisIO.wseis(InputDict["tmppath"]*"/"*fname_out, SremEQ)

    end

    return nothing
end

end
