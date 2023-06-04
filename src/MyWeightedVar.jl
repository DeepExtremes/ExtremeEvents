using WeightedOnlineStats: WeightedOnlineStat, smooth
import OnlineStatsBase # ._fit!, OnlineStatsBase._merge!
import Statistics
"""
    MyWeightedVariance(T = Float64)

Simple weighted variance, tracked as type `T`, ignoring values with 0 weight.

# Example:
    o = fit!(MyWeightedVariance(), rand(100), rand(100))
    sum(o)
    mean(o)
    var(o)
    std(o)
"""
mutable struct MyWeightedVariance{T} <: WeightedOnlineStat{T}
    μ::T
    σ2::T
    W::T
    W2::T
    n::Int
    function MyWeightedVariance{T}(
            μ = T(0), σ2 = T(0), W = T(0), W2 = T(0), n = 0
        ) where T
        new{T}(T(μ), T(σ2), T(W), T(W2), Int(n))
    end
end

MyWeightedVariance(μ::T, σ2::T, W::T, W2::T, n::Int) where T =
    MyWeightedVariance{T}(μ, σ2, W, W2, n)
MyWeightedVariance(::Type{T}) where T = MyWeightedVariance(T(0), T(0), T(0), T(0), 0)
MyWeightedVariance() = MyWeightedVariance(Float64)

function OnlineStatsBase._fit!(o::MyWeightedVariance{T}, x, w) where T
    xx = convert(T, x)
    ww = convert(T, w)

    # skip when ww == 0
    if ww == 0 return o end

    o.n += 1
    γ1 = T(1) / o.n
    o.W = smooth(o.W, ww, γ1)
    o.W2 = smooth(o.W2, ww * ww, γ1)
    # change γ, µ and σ2
    γ = ww / (o.W * o.n)
    μ = o.μ

    o.μ = smooth(o.μ, xx, γ)
    o.σ2 = smooth(o.σ2, (xx - o.μ) * (xx - μ), γ)

    return o
end

function OnlineStatsBase._merge!(
        o::MyWeightedVariance{T},
        o2::MyWeightedVariance
    ) where T

    o2_μ = convert(T, o2.μ)
    o2_σ2 = convert(T, o2.σ2)
    o2_W = convert(T, o2.W)
    o2_W2 = convert(T, o2.W2)


    n = o.n + o2.n
    W = smooth(o.W, o2_W, o2.n / n)
    γ1 = (o.W * o.n) / (W * n)
    γ2 = (o2_W * o2.n) / (W * n)

    μ = smooth(o.μ, o2_μ, γ2)

    # o.σ2 =
    #     γ1 * ( o.σ2 +  o.μ ^ 2) +
    #     γ2 * (o2_σ2 + o2_μ ^ 2) -
    #     μ ^ 2
    o.σ2 =
        γ1 * ( o.σ2  + (o.μ - μ) ^ 2) +
        γ2 * (o2_σ2 + (o2_μ - μ) ^ 2)

    o.n = n
    o.μ = μ
    o.W = W
    o.W2 = OnlineStats.smooth(o.W2, o2_W2, o2.n / o.n)


    ###########################################

    # μ = o.μ
    #
    # γ = o2_W / (o2_W + o.W)
    # δ = o2_μ - o.μ
    #
    # o.σ2 = OnlineStats.smooth(o.σ2, o2_σ2, γ) + (δ ^ 2) * γ * (1.0 - γ)
    # o.μ = OnlineStats.smooth(o.μ, o2_μ, γ)
    # # o.σ2 = o.σ2 + (o.W + o2_W) * (μ)
    #
    # o.W += o2_W
    # o.W2 += o2_W2

    return o
end

OnlineStatsBase.value(o::MyWeightedVariance) = o.σ2
Base.sum(o::MyWeightedVariance) = mean(o) * meanweight(o) * nobs(o)
Statistics.mean(o::MyWeightedVariance) = o.μ
function Statistics.var(
        o::MyWeightedVariance;
        corrected = true,
        weight_type = :analytic
    )
    if corrected
        if weight_type == :analytic
            value(o) / (1 - (o.W2 * nobs(o)) / (weightsum(o) ^ 2))
        elseif weight_type == :frequency
            value(o) / (weightsum(o) - 1) * weightsum(o)
        elseif weight_type == :probability
            error("If you need this, please make a PR or open an issue")
        else
            throw(ArgumentError("weight type $weight_type not implemented"))
        end
    else
        value(o)
    end
end
Statistics.std(o::MyWeightedVariance; kw...) = sqrt.(var(o; kw...))
Base.copy(o::MyWeightedVariance) = MyWeightedVariance(o.μ, o.σ2, o.W, o.W2, o.n)