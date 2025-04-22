# Sampling methods

using Random

"""
    cdf_sample(sampler)

Randomly select a `FlexLevel` from `sampler` using inverse transform sampling.

Also returns a "free" random number in [0, 1) for use in subsequent rejection sampling
(see [`rejection_sample`](@ref)).
"""
@inline function cdf_sample(sampler::FlexleSampler{Float64})
    if isempty(sampler.levels)
        throw(DomainError("no positive weights in FlexleSampler"))
    end
    local chosen_level::FlexLevel
    norm_rand_n = rand() * sampler.sum
    cum_sum = 0.0
    for level in sampler.levels
        cum_sum += level.sum
        chosen_level = level
        (cum_sum > norm_rand_n) && break
    end
    return chosen_level, (norm_rand_n - cum_sum + chosen_level.sum) / chosen_level.sum
end

@inline function cdf_sample(sampler::FlexleSampler{Int64})
    if isempty(sampler.levels)
        throw(DomainError("no positive weights in FlexleSampler"))
    end
    local chosen_level::FlexLevel
    norm_rand_n = rand(1:sampler.sum)
    cum_sum = 0
    for level in sampler.levels
        cum_sum += level.sum
        chosen_level = level
        (cum_sum >= norm_rand_n) && break
    end
    return chosen_level, rand()
end

"""
    rejection_sample(rand_n, level, weights)

Randomly select an index from `level` using rejection sampling given a starting `rand_n`
and a `Vector` of `weights`.

`rand_n` is generated in the course of inverse transform sampling (see [`cdf_sample`](@ref))
performed prior to rejection sampling.
"""
@inline function rejection_sample(rand_n::Float64, level::FlexLevel{T}, weights::Vector{T}) where {T<:WeightNumber}
    # l = Float64(length(level.elements))
    # while true
    #     r = rand_n * l
    #     i::Int64, r_floor::Float64 = fast_Int64(r)
    #     idx = level.elements[i + 1]   # +1 to offset for 1-indexing
    #     if weights[idx] > (r - r_floor) * level.max
    #         return idx
    #     end
    #     rand_n = rand()
    # end
    while true
        i = rand(level.elements)
        if weights[i] > rand_n * level.max
            return i
        end
        rand_n = rand()
    end
end
