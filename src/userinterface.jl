# User interface containing all exported Flexle functions

using Random

"""
    FlexleSampler{T}(weights)

Create a `FlexleSampler` from a `Vector` of `weights`.
"""
function FlexleSampler{T}(weights::Vector{T}) where {T<:WeightNumber}
    w_vector = Vector(weights)
    w_zero = zero(T)
    w_sum = w_zero

    if length(w_vector) == 0
        return FlexleSampler(Vector{FlexLevel{T}}(), w_vector, w_sum, Vector{Int64}(), nothing)
    end

    weights_nonzero = filter(x -> !iszero(x), weights)
    if isempty(weights_nonzero)
        return FlexleSampler(Vector{FlexLevel{T}}(), w_vector, w_sum, zeros(Int64, length(w_vector)), nothing)
    end
    
    w_min, w_max = typemax(T), w_zero
    for w in weights_nonzero
        (w < w_min) && (w_min = w)
        (w > w_max) && (w_max = w)
    end
    
    uppermost_log_bound = let 
        logmax = log2(w_max)
        ceil_logmax = ceil(logmax)
        (logmax == ceil_logmax) ? Int64(logmax) + 1 : Int64(ceil_logmax) # if w_max is a power of 2, can't take its ceiling to get its upper bound - need to add 1 instead
    end
    num_levels = uppermost_log_bound - floor_log2(w_min)     # e.g. -2,5 ==> 7 levels [4,3,2,1,0,-1,-2]
    
    levels = Vector{FlexLevel{T}}(undef, num_levels)   # add check for unreasonable number of levels before allocating space?
    element_positions = zeros(Int64, length(w_vector))

    lower_bound = lower_power_of_2_bound(w_min)
    for i in num_levels:-1:1
        upper_bound = lower_bound * 2
        levels[i] = FlexLevel{T}((lower_bound, upper_bound), w_zero, w_zero, 0, Vector{Int64}())
        lower_bound = upper_bound
    end

    for i in eachindex(w_vector)
        w::T = w_vector[i]
        iszero(w) && continue
        l = levels[level_index(w, uppermost_log_bound)]
        push!(l.elements, i)
        if w == l.max
            l.num_max += 1
        elseif w > l.max
            l.max = w
            l.num_max = 1
        end
        l.sum += w
        element_positions[i] = length(l.elements)
        w_sum += w
    end

    return FlexleSampler{T}(levels, w_vector, w_sum, element_positions, uppermost_log_bound)
end

function FlexleSampler(weights::Vector{T}) where {T<:WeightNumber}
    return FlexleSampler{T}(weights)
end

"""
    FlexleSampler()

Create an empty `FlexleSampler`.
"""
function FlexleSampler(; T=Float64)
    return FlexleSampler{T}(Vector{T}())
end

"""
    getindex(sampler, i)

Get the weight of element `i` in `sampler`.
"""
function Base.getindex(sampler::FlexleSampler{T}, i::Int64) where {T<:WeightNumber}
    return sampler.weights[i]
end

"""
    setindex!(sampler, w, i)

Set the weight of element `i` in `sampler` equal to `w`, returning the difference between the new and old values of `i`.
"""
function Base.setindex!(sampler::FlexleSampler{T}, w::T, i::Int64) where {T<:WeightNumber}
    from::Union{Nothing, FlexLevel{T}} = nothing
    to::Union{Nothing, FlexLevel{T}} = nothing
    levels = sampler.levels
    w_old::T = sampler.weights[i]
    delta::T = w - w_old
    nonzero = !iszero(w_old), !iszero(w)
    if nonzero[1]
        from = get_level(w_old, sampler)
    end
    if nonzero[2]
        bounds = logbounds(w)
        if !in_sampler(bounds, sampler)
            extend_levels!(bounds, sampler)
        end
        to = get_level(bounds, sampler)
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
        remove_from_FlexLevel!(i, from, sampler, update_sampler_sum=false)
    end
    sampler.weights[i] = w  # weight vector update must be between removal and addition
    if nonzero[2]
        add_to_FlexLevel!(i, to, sampler, update_sampler_sum=false)
    end
    sampler.sum += delta

    if nonzero[1] && (from === levels[begin] || from === levels[end]) && !level_is_populated(from)
        trim_trailing_levels!(sampler)
    end

    return delta
end

"""
    getweights(sampler)

Get the `Vector` of all the weights in `sampler`.
"""
function getweights(sampler::FlexleSampler{T}) where {T<:WeightNumber}
    return sampler.weights
end

"""
    numweights(sampler)

Get the number of weights in `sampler`.
"""
function numweights(sampler::FlexleSampler{T}) where {T<:WeightNumber}
    return length(sampler.weights)
end

"""
    push!(sampler, w)

Add a new element with weight `w` to `sampler`, updating all fields accordingly.

Returns the new number of weights in `sampler`, or equivalently, the index corresponding to the new element added.
"""
function Base.push!(sampler::FlexleSampler{T}, w::T) where {T<:WeightNumber}
    push!(sampler.weights, w)
    if !iszero(w)
        bounds = logbounds(w)
        if !in_sampler(bounds, sampler)
            extend_levels!(bounds, sampler)
        end
        to = get_level(bounds, sampler)
        add_to_FlexLevel!(length(sampler.weights), to, sampler)
    else
        push!(sampler.element_positions, 0)
    end
    return length(sampler.weights)
end

"""
    deleteat!(sampler, i)

Remove element `i` from `sampler` completely, updating all fields accordingly.

Returns the new number of weights in `sampler`.

All elements of index `>i` are updated to account for the removal of element `i`.
"""
function Base.deleteat!(sampler::FlexleSampler{T}, i::Int64) where {T<:WeightNumber}
    w::T = sampler.weights[i]
    if !iszero(w)
        bounds = logbounds(w)
        from = get_level(bounds, sampler)
        remove_from_FlexLevel!(i, from, sampler)
    end
    deleteat!(sampler.weights, i)
    for level in sampler.levels
        elements = level.elements
        for j in eachindex(elements)
            if elements[j] > i
                elements[j] -= 1
            end
        end
    end
    deleteat!(sampler.element_positions, i)
    if !iszero(w) && (from === sampler.levels[begin] || from === sampler.levels[end]) && isempty(from.elements)
        trim_trailing_levels!(sampler)
    end
    return length(sampler.weights)
end

"""
    sample(sampler)

Take a single random sample from `sampler`, returning the index of the element sampled.

Samples by inverse transform sampling (see `cdf_sample`(@ref)) to select a `FlexLevel` in `sampler`, then rejection
sampling (see `rejection_sample`(@ref)) an index from said `FlexLevel`.
"""
@inline function sample(sampler::FlexleSampler{T}) where {T<:WeightNumber}
    level, rand_n = cdf_sample(sampler)
    return rejection_sample(rand_n, level, sampler.weights)
end

