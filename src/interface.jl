# User interface containing all exported Flexle functions

using Revise
using Random

"""
    FlexleSampler(weights)

Create a `FlexleSampler` from a `Vector` of `weights`.
"""
function FlexleSampler(weights::AbstractVector{Float64})
    w_vector = Vector(weights)
    w_sum = 0.0

    if length(w_vector) == 0
        return FlexleSampler(Vector{FlexLevel}(), w_vector, w_sum, Vector{Int64}(), nothing)
    end

    weights_nonzero = filter(x -> !iszero(x), weights)
    if isempty(weights_nonzero)
        return FlexleSampler(Vector{FlexLevel}(), w_vector, w_sum, zeros(Int64, length(w_vector)), nothing)
    end
    
    w_min, w_max = Inf, 0.0
    for w in weights_nonzero
        (w < w_min) && (w_min = w)
        (w > w_max) && (w_max = w)
    end
    
    uppermost_log_bound = let 
        logmax = log2(w_max)
        ceil_logmax = ceil(logmax)
        (logmax == ceil_logmax) ? Int64(logmax) + 1 : Int64(ceil_logmax) # if w_max is a power of 2, can't take its ceiling to get its upper bound - need to add 1 instead
    end
    num_levels = uppermost_log_bound - floorLog2(w_min)     # e.g. -2,5 ==> 7 levels [4,3,2,1,0,-1,-2]
    
    levels = Vector{FlexLevel}(undef, num_levels)   # add check for unreasonable number of levels before allocating space?
    index_positions = zeros(Int64, length(w_vector))

    lower_bound = lowerPowerOf2Bound(w_min)
    for i in num_levels:-1:1
        upper_bound = lower_bound * 2.0
        levels[i] = FlexLevel((lower_bound, upper_bound), 0.0, 0.0, Vector{Int64}())
        lower_bound = upper_bound
    end

    for i in eachindex(w_vector)
        w::Float64 = w_vector[i]
        iszero(w) && continue
        l = levels[levelIndex(w, uppermost_log_bound)]
        push!(l.indices, i)
        (w > l.max) && (l.max = w)
        l.sum += w
        index_positions[i] = length(l.indices)
        w_sum += w
    end

    return FlexleSampler(levels, w_vector, w_sum, index_positions, uppermost_log_bound)
end

"""
    getindex(sampler, i)

Get the weight of element `i` in `sampler`.
"""
function Base.getindex(sampler::FlexleSampler, i::Int64)
    return sampler.weights[i]
end

"""
    setindex!(sampler, w, i)

Set the weight of element `i` in `sampler` equal to `w`, returning the difference between the new and old values of `i`.
"""
function Base.setindex!(sampler::FlexleSampler, w::Float64, i::Int64)
    from::Union{Nothing, FlexLevel} = nothing
    to::Union{Nothing, FlexLevel} = nothing
    levels = sampler.levels
    w_old::Float64 = sampler.weights[i]
    delta::Float64 = w - w_old
    nonzero = !iszero(w_old), !iszero(w)
    if nonzero[1]
        from = getLevel(w_old, sampler)
    end
    if nonzero[2]
        bounds = logBounds(w)
        if !inSampler(bounds, sampler)
            extendLevels!(bounds, sampler)
        end
        to = getLevel(bounds, sampler)
        if from == to
            sampler.weights[i] = w
            to.sum += delta
            sampler.sum += delta
            if w > to.max
                to.max = w
            end
            return delta
        end
    end

    if nonzero[1]
        removeFromFlexLevel!(i, from, sampler, update_sampler_sum=false)
    end
    sampler.weights[i] = w  # weight vector update must be between removal and addition
    if nonzero[2]
        addToFlexLevel!(i, to, sampler, update_sampler_sum=false)
    end
    sampler.sum += delta

    if nonzero[1] && (from === levels[begin] || from === levels[end]) && !levelIsPopulated(from)
        trimTrailingLevels!(sampler)
    end

    return delta
end

"""
    push!(sampler, w)

Add a new element with weight `w` to `sampler`, updating all fields accordingly.

Returns the new number of weights in `sampler`, or equivalently, the index corresponding to the new element added.
"""
function Base.push!(sampler::FlexleSampler, w::Float64)
    push!(sampler.weights, w)
    if !iszero(w)
        bounds = logBounds(w)
        if !inSampler(bounds, sampler)
            extendLevels!(bounds, sampler)
        end
        to = getLevel(bounds, sampler)
        addToFlexLevel!(length(sampler.weights), to, sampler)
    else
        push!(sampler.index_positions, 0)
    end
    return length(sampler.weights)
end

"""
    deleteat!(sampler, i)

Remove element `i` from `sampler` completely, updating all fields accordingly.

Returns the new number of weights in `sampler`.

All elements of index `>i` are updated to account for the removal of element `i`.
"""
function Base.deleteat!(sampler::FlexleSampler, i::Int64)
    w::Float64 = sampler.weights[i]
    if !iszero(w)
        bounds = logBounds(w)
        from = getLevel(bounds, sampler)
        removeFromFlexLevel!(i, from, sampler)
    end
    deleteat!(sampler.weights, i)
    for level in sampler.levels
        indices = level.indices
        for j in eachindex(indices)
            if indices[j] > i
                indices[j] -= 1
            end
        end
    end
    deleteat!(sampler.index_positions, i)
    if !iszero(w) && (from === sampler.levels[begin] || from === sampler.levels[end]) && isempty(from.indices)
        trimTrailingLevels!(sampler)
    end
    return length(sampler.weights)
end

"""
    sample(sampler)

Take a single random sample from `sampler`, returning the index of the element sampled.

Samples by inverse transform sampling (see `cdfSample`(@ref)) to select a `FlexLevel` in `sampler`, then rejection
sampling (see `rejectionSample`(@ref)) an index from said `FlexLevel`.
"""
@inline function sample(sampler::FlexleSampler)
    level, rand_n = cdfSample(sampler)
    return rejectionSample(rand_n, level, sampler.weights)
end

