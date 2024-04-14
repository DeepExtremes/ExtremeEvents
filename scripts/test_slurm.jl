
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

y = 1950; p = [y:(y+12)]

while y < 2010
    global y += 10
    append!(p,[y:(y+12)])
end
# period = p
# aperiod = map( x -> replace(string(x), ":" => "_", r"^" => "_"),p)
@everywhere begin
    land = ""
    region = ""
end

@sync @distributed for period in p
    @show aperiod = replace(string(period), ":" => "_", r"^" => "_")
    outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/test_$aperiod$land$region.txt"
    txt = "This is $period"
    open(outpath, "w") do f
        write(f,txt)
    end
end