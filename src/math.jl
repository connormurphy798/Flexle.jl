const EXPONENT_MASK_FLOAT64::Int64      = 0x7FF0000000000000
const EXPONENT_SHIFT_FLOAT64::Int64     = 52
const EXPONENT_OFFSET_FLOAT64::Int64    = 1023
const MANTISSA_MASK_FLOAT64::Int64      = 0x000FFFFFFFFFFFFF
const MANTISSA_PLUS1_FLOAT64::Int64     = 0x0010000000000000

import Base.&

Base.:&(f::Float64, i::Int64) = reinterpret(Int64, f) & i

"""
    get_exponent_bits(a)

Get an Int64 representing the exponent bits in the IEEE754 representation of `a`.
"""
@inline function get_exponent_bits(a::Float64)
    return (a & EXPONENT_MASK_FLOAT64) >> EXPONENT_SHIFT_FLOAT64
end

"""
    fast_Int64(n)

Truncate and convert `n` to an `Int64`, but only safely for positive numbers.

Returns the `Int64`-converted and truncated values of `n` as a 2-tuple.
"""
@inline function fast_Int64(n::Float64)
    exponent_bits = get_exponent_bits(n) - EXPONENT_OFFSET_FLOAT64
    shift = EXPONENT_SHIFT_FLOAT64 - exponent_bits
    return ((n & MANTISSA_MASK_FLOAT64) | MANTISSA_PLUS1_FLOAT64) >> shift, reinterpret(Float64, (reinterpret(Int64, n) >> shift) << shift)
end

"""
    approxeq(a, b, t=1e-9)

Return a `Bool` indicating whether `a` and `b` are equal to within `t`.
"""
function approxeq(a::Float64, b::Float64; t::Float64=1e-9)
    return abs(a-b) < t
end

"""
    lower_power_of_2_bound(n)

Get the largest power of 2 less than or equal to a non-negative `n`.

# Examples

`lower_power_of_2_bound(33.0)` ==> `32.0`

`lower_power_of_2_bound(32.0)` ==> `32.0`

`lower_power_of_2_bound(0.75)` ==> `0.5`
"""
function lower_power_of_2_bound(n::Float64)
    return reinterpret(Float64, n & EXPONENT_MASK_FLOAT64)
end

"""
    floor_log2(n)

Get the `Int64` floor of the log2 of `n`.

# Examples

`floor_log2(33.0)` ==> `5`

`floor_log2(32.0)` ==> `5`

`floor_log2(0.75)` ==> `-1`
"""
function floor_log2(n::Float64)
    return ((n & EXPONENT_MASK_FLOAT64) >> EXPONENT_SHIFT_FLOAT64) - EXPONENT_OFFSET_FLOAT64
end

"""
    logbounds(n)

Return a tuple `l,u` giving two adjacent powers of 2 such that `l <= n < u`.
"""
function logbounds(n::Float64)
    l = lower_power_of_2_bound(n)
    return l, l*2.0
end

"""
    level_index(w, u)

Given a weight `w`, return the index of the level in some `FlexleSampler.levels` with maximum upper bound `2^u` where `w` would belong.

Returns `0` if `w` is `0.0`, indicating that `w` belongs in no level.

# Examples

`level_index(14.2, 6)` ==>  `3`

The value `14.2` belongs in the `(8.0, 16.0)` level, which in a `FlexleSampler` that
starts with a level of bounds `(32.0, 64.0)` (`64` being `2^6`) is at `levels[3]`.

`level_index(8.0, 4)` ==> `1`

`level_index(8.0, 5)` ==> `2`
"""
function level_index(w::Float64, u::Union{Int64, Nothing})
    return iszero(w) || isnothing(u) ? 0 : u - floor_log2(w)
end

"""
    level_index(bounds, u)

Get the index of the `FlexLevel` that has bounds given by `bounds` in some `FlexleSampler.levels` with maximum upper bound `2^u`.

Returns `0` if no such level exists.
"""
function level_index(bounds::Tuple{Float64,Float64}, u::Union{Int64, Nothing})
    return level_index(bounds[1], u)
end

"""
    get_level(bounds, sampler)

Return the `FlexLevel` in `sampler.levels` with bounds given by `bounds`.

Returns `nothing` if no such level exists.
"""
function get_level(bounds::Tuple{Float64,Float64}, sampler::FlexleSampler)
    l = level_index(bounds[1], sampler.max_log2_upper_bound)
    return l==0 ? nothing : sampler.levels[l]
end

"""
    get_level(w, levels)

Return the `FlexLevel` in `sampler.levels` where a weight `w` would belong.

Returns `nothing` if no such level exists.
"""
function get_level(w::Float64, sampler::FlexleSampler)
    l = level_index(w, sampler.max_log2_upper_bound)
    return l==0 ? nothing : sampler.levels[l]
end

"""
    log_dist(a, b)

Return the floor of the log2 of `a/b`. 
"""
function log_dist(a::Float64, b::Float64)
    return floor_log2(b) - floor_log2(a)
end

"""
    in_sampler(i, sampler)

Return a `Bool` indicating whether an index `i` is present at some `FlexLevel` in `sampler`.
"""
function in_sampler(i::Int64, sampler::FlexleSampler)
    for level in sampler.levels
        if i in level.elements
            return true
        end
    end
    return false
end

"""
    in_sampler(bounds, sampler)

Return a `Bool` indicating whether a `FlexLevel` with bounds `bounds` is present in `sampler`.
"""
function in_sampler(bounds::Tuple{Float64,Float64}, sampler::FlexleSampler)
    return (length(sampler.levels) > 0 && sampler.levels[begin].bounds[1] >= bounds[1]) && (bounds[1] >= sampler.levels[end].bounds[1]) # bounds between largest and smallest levels' bounds (inclusive)
end

"""
    level_is_populated(level)

Return a `Bool` indicating whether `level` contains any elements.
"""
function level_is_populated(level::FlexLevel)
    return !isempty(level.elements)
end

"""
    max_level_weight(level, sampler)

Get the largest weight in `sampler` of any element in `level` and the number of times it occurs.
"""
function max_level_weight(level::FlexLevel, sampler::FlexleSampler)
    m::Float64 = 0.0
    weights::Vector{Float64} = sampler.weights
    n::Int64 = 0
    for i in level.elements
        w = weights[i]
        if w == m
            n += 1
        elseif w > m
            m = w
            n = 1
        end
    end
    return m, n
end

"""
    recalculate_sum
"""
function recalculate_sum(w::Float64, sum::Float64; tol=0.999) # if one weight is 99.9% of sum, all other weights have >=10 digits of lost precision
    recalc = !iszero(sum) && w/sum >= tol                      # log2( 99.9 / (100-99.9) ) = 9.96
    return recalc
end
