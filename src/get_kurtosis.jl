module Get_kurtosis

export get_kurtosis

using Distributions, SeisIO

"""
    get_kurtosis(data::SeisData,timewinlength::Float64=60)

    compute kurtosis at each timewindow

# Input:
    - `data::SeisData`    : SeisData from SeisIO
    - `timewinlength::Float64`  : time window to calculate kurtosis

    kurtosis evaluation following Baillard et al.(2013)
"""

function get_kurtosis(data::SeisChannel, timewinlength::Float64=60)

    #convert window lengths from seconds to samples
    TimeWin = trunc(Int,timewinlength * data.fs)

    #set long window length to user input since last window of previous channel will have been adjusted
    TimeWin = trunc(Int,timewinlength * data.fs)
    data.misc["kurtosis"] = zeros(Float64, length(data.x))

    #reset long window counter and triggers for current channel
    trigger = Float64[]
    i = 0

    #loop through current channel by sliding
    while i < length(data.x) - TimeWin
        #define chunk of data based on long window length and calculate long-term average
        Trace = @views data.x[i+1:i+TimeWin]
        data.misc["kurtosis"][i+TimeWin] = kurtosis(Trace)
        #advance time window
        i += 1
    end

    return data

end

end
