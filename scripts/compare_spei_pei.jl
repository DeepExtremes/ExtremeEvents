#This reads pei dataand computes daily SPEI based on fitting a log-logistic distribution
#Compares thresholds derived from quantiles of the data as well as the estimated distribution

using YAXArrays
using CairoMakie
using SpecialFunctions: gamma
import NetCDF
using NonlinearSolve: NonlinearProblem, solve, NewtonRaphson, SimpleNewtonRaphson, IntervalNonlinearProblem
using Statistics: quantile, quantile!, Statistics
using Distributions: Normal
ds = open_dataset("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/PEICube.zarr")

g(fi,di,s) = (1-fi)^s * di
function w(s,di) 
    N = length(di)
    agg = 0.0
    for i in 1:N
        agg = agg+g((i-0.35)/N,di[i],s)
    end
    agg/N
end
function estimate_params(data)
    di = sort!(data)    
    w0 = w(0,di)
    w1 = w(1,di)
    w2 = w(2,di)
    β = (2*w1 - w0) / (6*w1 - w0 - 6*w2)
    α = (w0 - 2*w1)*β / gamma(1 + 1/β) / gamma(1-1/β)
    γ = w0 - α * gamma(1 + 1/β) * gamma(1-1/β)
    α, β, γ
end

struct SPEI
    α::Float64
    β::Float64
    γ::Float64
end
function (spei::SPEI)(x)
    x < spei.γ && return NaN
    spei.β < 0.0 && return NaN
    spei.α < 0.0 && return NaN
    any(isnan,(spei.α, spei.β, spei.γ)) && return NaN
    F = inv(1 + (spei.α / (x-spei.γ))^spei.β)
    P = 1-F
    C0 = 2.515517; C1 = 0.802853; C2 = 0.010328; d1 = 1.432788; d2 = 0.189269; d3 = 0.001308
    if P<=0.5
        W = sqrt(-2*log(P))
        return W - (C0 + C1*W + C2*W*W)/(1 + d1*W + d2*W*W + d3*W*W*W)
    else
        W = sqrt(-2*log(1-P))
        return -W  + (C0 + C1*W + C2*W*W)/(1 + d1*W + d2*W*W + d3*W*W*W)
    end
end
pdf(x,sp::SPEI) = sp.β / sp.α * ((x-sp.γ)/sp.α)^(sp.β-1) * (1 + ((x-sp.γ)/sp.α)^sp.β)^(-2)
cdf(x,sp::SPEI) = inv(1 + ((x-sp.γ)/sp.α)^(-sp.β))
Statistics.quantile(sp::SPEI,q::Number) = sp.α * (q/(1-q))^(1/sp.β) + sp.γ
Statistics.quantile(sp::SPEI,q) = map(Base.Fix1(quantile,sp),q)
speifunc(x,(spei,q)) = spei(x)-q
speifunc2(x,(spei,q)) = cdf(x,spei)-q
isgooddata(x) = !ismissing(x) && !isnan(x)
function normalize_pei(xout,parout,qout,xin)
    gooddata = findall(isgooddata,xin)
    xinclean = collect(Float64,filter(isgooddata,xin))
    sort!(xinclean)
    qs = (0.01,0.05,0.1) 
    qout[2,:] .= quantile!(xinclean,qs)
    α, β, γ = estimate_params(xinclean)
    parout[1] = α
    parout[2] = β
    parout[3] = γ
    spei = SPEI(α,β,γ)
    for i in 1:3
        q = qs[i]
        qout[1,i] = quantile(spei,q)
    end
    fill!(xout,NaN)
    for i in gooddata
        xout[i] = spei(xin[i])
    end
end
dssub = ds.pei_30
quantdim = YAXArrays.Dim{:quantile}([0.01,0.05,0.1])
quantmethod = YAXArrays.Dim{:method}(["Distribution","Empirical"])
pardim = YAXArrays.Dim{:param}(["α","β","γ","success"])

r = mapCube(normalize_pei,dssub,
    indims=InDims("Ti"),
    outdims=(OutDims("Ti", path="/Net/Groups/BGI/work_2/scratch/fgans/SPEI_30.zarr",overwrite=true),
    OutDims(pardim,path="/Net/Groups/BGI/work_2/scratch/fgans/SPEI_30_params.zarr",overwrite=true),
    OutDims(quantmethod,quantdim,path="/Net/Groups/BGI/work_2/scratch/fgans/pei_30_quantiles.zarr",overwrite=true)),
    )


