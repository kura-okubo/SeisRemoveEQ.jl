#using Distributed
#addprocs(3)

@everywhere using SeisRemoveEQ, SeisIO

using Dates, JLD2

#==================================================#

# Use SeisDownload.jl to make jld2 file.
finame = "./dataset/3days_CInetwork_ridgecrest_fortest.jld2"

MAX_MEM_PER_CPU = 2.0 # [GB] maximum allocated memory for one cpu

#---Kurtosis parameters---#
IsKurtosisRemoval = true

max_edgetaper_duration = 60 * 5 #[s] tapering to avoid edge effect of response removal
kurtosis_timewindow = 60 * 3 #300 # timewindow to compute kurtosis [s]: large number causes longer computational time.
kurtosis_threshold = 3.0 #kurtosis threshold to detect outlier

#---STA/LTA parameters---#
IsSTALTARemoval = true

stalta_longtimewindow = 60*60*24 #[s] long time window for STA/LTA: Note that short time window is same with removal_shorttimewindow
stalta_threshold = 0.8 #STA/LTA threshold to detect earthquake and tremors: For our purpose, this threshold is smaller than ordinal detection

#---EQ removal parameter---#
invert_tukey_α = 0.01 #invert_tukey window (0<α<1: 1.0 for smooth removal but causes incomplete tremors removal)
removal_shorttimewindow = 60 * 3 # timewindow to remove EQ: removing data if this tw contains data more than threshold.
overlap = 60; # overlap of removal window[s]

#---Plotting variables---#
IsSaveFig = true # save figure to check removal
plot_kurtosis_α = 1.2 # amplitude normalization for plotting kurtosis
plot_boxheight = 1.5
plot_span = 100 # to reduce number of plot point

#---Output file info---#
fodir       = "./dataset"
foname      = "3days_CInetwork_ridgecrest_remEQ" # data is saved at ./dataset/$foname.jld2

#==================================================#

#intended to parallelization
InputDictionary = Dict( "MAX_MEM_PER_CPU"       => float(MAX_MEM_PER_CPU),
                        "finame"                => finame,
                        "IsKurtosisRemoval"     => IsKurtosisRemoval,
                        "max_edgetaper_duration" => float(max_edgetaper_duration),
                        "kurtosis_timewindow"   => float(kurtosis_timewindow),
                        "kurtosis_threshold"    => float(kurtosis_threshold),
                        "IsSTALTARemoval"       => IsSTALTARemoval,
                        "stalta_longtimewindow" => float(stalta_longtimewindow),
                        "stalta_threshold"      => float(stalta_threshold),
                        "invert_tukey_α"        => float(invert_tukey_α),
                        "removal_shorttimewindow"=> float(removal_shorttimewindow),
                        "overlap"               => float(overlap),
                        "plot_kurtosis_α"       => float(plot_kurtosis_α),
                        "plot_boxheight"        => float(plot_boxheight),
                        "plot_span"             => trunc(Int, plot_span),
                        "fodir"                 => fodir,
                        "foname"                => foname,
                        "IsSaveFig"             => IsSaveFig)

seisremoveEQ(InputDictionary)
