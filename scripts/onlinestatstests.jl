tmp = stack(DataFrame(a = rand(0:16, 20), b = rand(0:16, 20), cont = rand(1:3,20), lat = 80:-4:1, lsm=rand(20)), [:a, :b], value_name = :event)

function fit1(df)
    dfg = groupby(df, [:variable, :cont])
    continents = range(1,3)
    allhists = Dict("$(i).$(k)" => WeightedHist(-0.5:1.0:16.5) for i in ["a","b"], k = continents)
    for k in keys(dfg)
        if !ismissing(k[2])
            dfs = dfg[k]
            fit!(allhists["$(k[1]).$(k[2])"], dfs.event, cosd.(dfs.lat) .* (dfs.lsm .> 0.5))
        end
    end
    allhists
end

mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))

tst = fit1(tmp)

t = (tmp[1:10,:], tmp[11:20,:], tmp[21:30,:], tmp[31:40,:])
tst2 = pmapreduce(mergefun,t) do tab
    fit1(tab)
end