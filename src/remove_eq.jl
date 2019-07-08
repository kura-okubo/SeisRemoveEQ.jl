module Remove_eq

export detect_eq, remove_eq

using DSP, Distributions, Statistics, SeisIO, Printf, PlotlyJS


"""
    detect_eq(data::SeisData,tw::Float64=60.0, threshold::Float64=3.0, overlap::Float64=30)

find earthquake by kurtosis threshold

# Input:
    - `data::SeisData`    : SeisData from SeisIO
    - `tw::Float64`  : time window to evaluate if earthquake is contained.
    - `threshold::Float64` : kurtosis threshold: if kurtosis > threshold, the time window contains earthquake
    - `overlap::Float64`           : overlap of time window to control the margin of earthquake removal. (large overlap assigns large margin before and after earthquake.)

    kurtosis evaluation following Baillard et al.(2013)
"""
function detect_eq(data::SeisChannel,tw::Float64=60.0, kurtosis_threshold::Float64=3.0, overlap::Float64=30)

    #convert window lengths from seconds to samples
    twsize = trunc(Int, tw * data.fs)
    overlapsize = trunc(Int, overlap * data.fs)

    #calculate how much to move beginning of window each iteration
    slide = twsize-overlapsize

    #set long window length to user input since last window of previous channel will have been adjusted
    data.misc["eqtimewindow"] = fill(false, length(data.x))

    #kurtosis of timeseries
    ku1 = data.misc["kurtosis"][:]
    #reset long window counter and triggers for current channel
    i = 0

    #loop through current channel by sliding
    while i < length(ku1) - twsize

        #check if last window and change long window length if so
        if length(ku1) - i < twsize
            twsize = length(ku1)-i
        end

        #define chunk of data based on long window length and calculate long-term average
        twtrace = @views ku1[i+1:i+twsize]

        if !isnothing(findfirst(x -> x > kurtosis_threshold, twtrace))
            #this time window includes earthquake
            for tt= i+1:i+twsize
                data.misc["eqtimewindow"][tt] = true
            end
        end

        #advance long window
        i = i + slide

    end


    return data

end


"""
    remove_eq(data::SeisData)

remove earthquake by kurtosis and STA/LTA threshold

# Input:
    - `data::SeisData`    : SeisData from SeisIO

"""
function remove_eq(data::SeisChannel, data_origin::SeisChannel, invert_tukey_α::Float64, plot_kurtosis_α::Float64,
    plot_boxheight::Float64, plot_span::Int64, fodir::String, tstamp::String, tvec::Array{Float64,1}, IsSaveFig::Bool)

    eqidlist = data.misc["eqtimewindow"][:]

    i = 1

    t1 = []
    t2 = []
    y1 = []
    y2 = []

    while i <= length(data.x)
        if eqidlist[i]
            push!(t1, tvec[i])

            t1id = i
            #find next id
            t2id = t1id + findfirst(x -> x == false, eqidlist[t1id:end])

            if t2id == t1id; t2id = length(eqidlist); end;

            push!(t2, tvec[t2id])

            # apply invert tukey window
            invtukeywin = -tukey(t2id-t1id+1, invert_tukey_α) .+ 1

            data.x[t1id:t2id] = data.x[t1id:t2id] .* invtukeywin

            iinc = t2id - i

            #boxsize
            push!(y1, -plot_boxheight)
            push!(y2, plot_boxheight)

        else
            iinc = 1
        end
        i += iinc
    end

    if IsSaveFig

        normalized_amp = 0.5 * maximum(data_origin.x[1:plot_span:end])

        trace1 = scatter(;x=tvec[1:plot_span:end], y=data_origin.x[1:plot_span:end]./ normalized_amp,
         mode="lines", name="raw data", line=attr(line_color="black", line_width=1))

        trace2 = scatter(;x=tvec[1:plot_span:end], y=data.x[1:plot_span:end]./ normalized_amp,
         mode="lines", name="after remove", line=attr(line_color="blue", line_width=2))

        trace3 = scatter(;x=tvec[1:plot_span:end], y=data.misc["kurtosis"][1:plot_span:end] ./ maximum(abs.(data.misc["kurtosis"][1:plot_span:end])) .* plot_kurtosis_α,
         line = attr(color="red", line_width=1), mode="lines", name="kurtosis")

        if !isempty(t1)
         shapes = PlotlyJS.rect(t1, t2, y1, y2; fillcolor="#ff99cc", opacity=0.3, line_width=0)
         layout = Layout(shapes=shapes, width=1200, height=600,
             xaxis=attr(title="Time [hour]"),
             yaxis=attr(title="Normalized velocity"),
             font =attr(size=13),
             showlegend=true,
             title = @sprintf("%s %s", data.id, tstamp))

         p = plot([trace1; trace2; trace3],layout)
        else
         layout = Layout(width=1200, height=600,
             xaxis=attr(title="Time [hour]"),
             yaxis=attr(title="Normalized velocity"),
             font =attr(size=12),
             showlegend=true,
             title = @sprintf("%s %s", data.id, tstamp))

         p = plot([trace1; trace2; trace3],layout)

        end

        figdir = joinpath(fodir, "fig")
        mkpath(figdir)
        figname = @sprintf("%s/%s_%s.png", figdir, data.id, tstamp)
        savefig(p, figname)
        #display(p)
        #println("press return for next plot...")
        #readline()
    end

    return data

end


end
