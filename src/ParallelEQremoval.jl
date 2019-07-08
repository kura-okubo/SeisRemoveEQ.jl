module EQremoval

export ParallelEQremoval

using Dates, JLD2, SeisIO, DSP, PlotlyJS, Printf, ORCA, BenchmarkTools

include("get_kurtosis.jl")
include("remove_eq.jl")
using .Get_kurtosis, .Remove_eq

"""
    ParallelEQremoval(dlid, InputDict::Dict)
    remove earthquake and save it into jld2 file.

"""
function ParallelEQremoval(dlid, InputDict::Dict)

    #store data

    finame               =InputDict["finame"]
    maxedgetaperduration =InputDict["maxedgetaperduration"]
    kurtosis_timewindow  =InputDict["kurtosis_timewindow"]
    threshold            =InputDict["threshold"]
    overlap              =InputDict["overlap"]
    tukey_α              =InputDict["tukey_α"]
    kurtosis_α           =InputDict["kurtosis_α"]
    boxheight            =InputDict["boxheight"]
    span                 =InputDict["span"]
    boxheight            =InputDict["boxheight"]
    IsSaveFig            =InputDict["IsSaveFig"]


    t = jldopen(finame)
    DLtimestamplist = t["info/DLtimestamplist"]
    stationlist = t["info/stationlist"]

    tstamp = DLtimestamplist[dlid]

    if mod(dlid, round(0.1*length(DLtimestamplist))+1) == 0
        println("start timestamp id: $dlid")
    end

    Eall = SeisData(length(stationlist))

    icount = 0

    bt_getkurtosis = 0.0
    bt_removeeq = 0.0

    for st = stationlist
        S = t[joinpath(tstamp, st)]

        if S.misc["dlerror"] == 0
            dt = 1/S.fs
            tvec = collect(0:S.t[2,1]-1) * dt ./ 60 ./ 60

            #tapering to avoid instrumental edge artifacts
            SeisIO.taper!(S,  t_max = maxedgetaperduration, α=0.05)
            S1 = deepcopy(S)


            bt_1 = @elapsed F1 = Get_kurtosis.get_kurtosis(S1, float(kurtosis_timewindow))
            tw = kurtosis_timewindow; #so far this is most stable
            bt_2 = @elapsed E1 = Remove_eq.remove_eq(F1, float(tw), float(threshold), float(overlap))

            bt_getkurtosis += bt_1
            bt_removeeq += bt_2

            eqidlist = E1[1].misc["eqtimewindow"][:]


            i = 1

            t1 = []
            t2 = []
            y1 = []
            y2 = []

            while i <= E1[1].t[2,1]
                if eqidlist[i]
                    push!(t1, tvec[i])

                    t1id = i
                    #find next id
                    t2id = t1id + findfirst(x -> x == false, eqidlist[t1id:end])

                    if t2id == t1id; t2id = length(eqidlist); end;

                    push!(t2, tvec[t2id])

                    # apply invert tukey window
                    invtukeywin = -tukey(t2id-t1id+1, tukey_α) .+ 1

                    E1[1].x[t1id:t2id] = E1[1].x[t1id:t2id] .* invtukeywin

                    iinc = t2id - i

                    #boxsize
                    push!(y1, -boxheight)
                    push!(y2, boxheight)

                else
                    iinc = 1
                end
                i += iinc
            end

            if IsSaveFig

                normalized_amp = 0.5 * maximum(S.x[1:span:end])

                trace1 = scatter(;x=tvec[1:span:end], y=S.x[1:span:end]./ normalized_amp,
                 mode="lines", name="raw data", line=attr(line_color="black", line_width=1))

                trace2 = scatter(;x=tvec[1:span:end], y=E1[1].x[1:span:end]./ normalized_amp,
                 mode="lines", name="after remove", line=attr(line_color="blue", line_width=2))

                trace3 = scatter(;x=tvec[1:span:end], y=F1[1].misc["kurtosis"][1:span:end] ./ maximum(abs.(F1[1].misc["kurtosis"][1:span:end])) .* kurtosis_α,
                 line = attr(color="red", line_width=1), mode="lines", name="kurtosis")

                if !isempty(t1)
                 shapes = PlotlyJS.rect(t1, t2, y1, y2; fillcolor="#ff99cc", opacity=0.3, line_width=0)
                 layout = Layout(shapes=shapes, width=1200, height=600,
                     xaxis=attr(title="Time [hour]"),
                     yaxis=attr(title="Normalized velocity"),
                     font =attr(size=13),
                     showlegend=true,
                     title = @sprintf("%s %s", tstamp, st))

                 p = plot([trace1; trace2; trace3],layout)
                else
                 layout = Layout(width=1200, height=600,
                     xaxis=attr(title="Time [hour]"),
                     yaxis=attr(title="Normalized velocity"),
                     font =attr(size=12),
                     showlegend=true,
                     title = @sprintf("%s %s", tstamp, st))

                 p = plot([trace1; trace2; trace3],layout)

                end

                mkpath("./fig")
                figname = @sprintf("./fig/%s_%s.png", tstamp, st)
                savefig(p, figname)
                #display(p)
                #println("press return for next plot...")
                #readline()
            end

        else
            #download error found: save as it is.
            E1 = SeisData(S)
        end
        icount += 1
        Eall[icount] = E1[1]
    end

    return (Eall, bt_getkurtosis, bt_removeeq)
end

end
