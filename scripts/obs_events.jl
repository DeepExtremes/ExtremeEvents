using DataFrames, Dates #DateFormats
import CSV

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/v3"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

# load list of events
obs0 = CSV.read("$path/../EventPart1_csv_Ltime.csv", DataFrame)
rename!(obs0, Dict("Name" => :Area, 
    "Event type" => :Event, 
    "when_from" => :Start, 
    "when_until" => :End, 
    "where_SW" => :West, 
    "where_SE" => :East,
    "where_NW" => :South,
    "where_NE" => :North,)
)
obs0.Start .= replace.(obs0.Start, r"\." => "-");
obs0.End .= replace.(obs0.End, r"\." => "-");
# remove spaces in columns names. DataConvenience.cleannames! from DataConvenience.jl
sort!(obs0, :Start, rev=true)
# select drought and heatwave (drop floods)
filter!(:Event=>!=("flood"), obs0)
# drop Continent
select!(obs0, Not(:Continent, :Year))

### EventPart2
obs = CSV.read("$(path)../EventPart3.csv", DataFrame; header=3)
# obs = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart2.csv", DataFrame; header=3)
# clean obs
obs = dropmissing!(obs)
obs.Start .= replace.(obs.Start, r"\." => "-");
obs.End .= replace.(obs.End, r"\." => "-");

# add obs_event
obs0.obs_event .= nrow(obs) .+ (1:nrow(obs0))
# reorder columns
select!(obs0, :obs_event, :)


# merge 2 tables:
obs = vcat(obs, obs0)

# sort by starting date
sort(obs, :Start)

# export to latex table
show(stdout, MIME("text/latex"),select(obs, Not(:obs_event)))
CSV.write(path * "refEvents.csv", obs)
