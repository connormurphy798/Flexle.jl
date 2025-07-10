# Structs and maintenance functions internal to the module

using Printf

# Maintenance methods

"""
    add_to_FlexLevel!(i, level, sampler, update_sampler_sum=true)

Place the element of index `i` in `sampler` into `level`.

Also updates `level.sum` and (iff `update_sampler_sum`) `sampler.sum` to reflect the addition of the weight of element `i`
to `level` and `sampler`, respectively.

# Note on `update_sampler_sum`

Keyword argument `update_sampler_sum` should only be `false` when `add_to_FlexLevel!` is used in conjunction with
[`remove_from_FlexLevel!`](@ref) to update the weight of an existing element `i` in `sampler`, as in the
following call pattern: 
```
remove_from_FlexLevel!(i, old_level, sampler, update_sampler_sum=false)
add_to_FlexLevel!(i, new_level, sampler, update_sampler_sum=false)
sampler.sum += new_i_weight - old_i_weight
``` 
This option is provided for performance purposes, as reading and writing to `sampler.sum` can be expensive. The caller MUST update
`sampler.sum` themselves if calling with `update_sampler_sum=false`.
"""
function add_to_FlexLevel!(i::Int64, level::FlexLevel, sampler::FlexleSampler; update_sampler_sum::Bool=true)
    push!(level.elements, i)
    w::Float64 = sampler.weights[i]
    level.sum += w
    (update_sampler_sum) && (sampler.sum += w)
    if w > level.max
        level.max = w
        level.num_max = 1
        
    elseif w == level.max
        level.num_max += 1
    end
    if i <= length(sampler.element_positions)
        sampler.element_positions[i] = length(level.elements)
    else
        push!(sampler.element_positions, length(level.elements))
    end
end

"""
    remove_from_FlexLevel!(i, level, sampler, update_sampler_sum=true)

Remove the element of index `i` in `sampler` from `level`.

Also updates `level.sum` and (iff `update_sampler_sum`) `sampler.sum` to reflect the removal of the weight of element `i`
from `level` and `sampler`, respectively.

# Note on `update_sampler_sum`

Keyword argument `update_sampler_sum` should only be `false` when `remove_from_FlexLevel!` is used in conjunction with
[`add_to_FlexLevel!`](@ref) to update the weight of an existing element `i` in `sampler`, as in the following call pattern: 
```
remove_from_FlexLevel!(i, old_level, sampler, update_sampler_sum=false)
add_to_FlexLevel!(i, new_level, sampler, update_sampler_sum=false)
sampler.sum += new_i_weight - old_i_weight
``` 
This option is provided for performance purposes, as reading and writing to `sampler.sum` can be expensive. The caller MUST update
`sampler.sum` themselves if calling with `update_sampler_sum=false`.
"""
function remove_from_FlexLevel!(i::Int64, level::FlexLevel, sampler::FlexleSampler; update_sampler_sum::Bool=true)
    w::Float64 = sampler.weights[i]
    len = length(level.elements)
    idx = sampler.element_positions[i]
    last = pop!(level.elements)
    if idx != len   # take last index and put it in the place of the one to be removed, unless the last one is itself to be removed
        level.elements[idx] = last
        sampler.element_positions[last] = idx
    end
    sampler.element_positions[i] = 0
    level.sum -= w
    (update_sampler_sum) && (sampler.sum -= w)
    if !level_is_populated(level)
        level.max = 0.0
        level.num_max = 0
    elseif w == level.max
        level.num_max -= 1
        if level.num_max == 0
            level.max, level.num_max = max_level_weight(level, sampler)
        end
    end
end

"""
    update_within_same_FlexLevel!(i, w, delta, level, sampler)

Update the weight `w` of element `i` within the appropriate `level` of `sampler`.

Only to be called when the prior weight of `i` in `sampler.weights[i]` belongs in the same `level` as new weight `w`.
"""
function update_within_same_FlexLevel!(i::Int64, w::Float64, delta::Float64, level::FlexLevel, sampler::FlexleSampler)
    w_old = sampler.weights[i]
    sampler.weights[i] = w
    level.sum += delta
    sampler.sum += delta
    if w_old == level.max
        level.num_max -= 1
        if level.num_max == 0
            level.max, level.num_max = max_level_weight(level, sampler)
        end
    elseif w > level.max
        level.max = w
        level.num_max = 1
    end
end

"""
    extend_levels!(bounds, sampler)

Extend `sampler.levels` to contain all appropriate `FlexLevel`s up to and including that specified by `bounds`.

Throws an error if a level with such bounds already exists in `levels`.
"""
function extend_levels!(bounds::Tuple{Float64,Float64}, sampler::FlexleSampler)
    if bounds[1] * 2.0 != bounds[2]
        throw("Invalid bounds - must be two adjacent powers of 2")
    end

    l_bound = bounds[1]
    if length(sampler.levels) == 0
        num_new_levels = 1
        post = Vector{FlexLevel}(undef, num_new_levels)
        u_bound = l_bound * 2.0
        post[1] = FlexLevel((l_bound, u_bound), 0.0, 0.0, 0, Vector{Int64}())
        l_bound = u_bound
        append!(sampler.levels, post)
    else
        extend_up = l_bound > sampler.levels[begin].bounds[1]
        extend_down = l_bound < sampler.levels[end].bounds[1]
        if extend_up
            num_new_levels = log_dist(sampler.levels[begin].bounds[1], l_bound)
            pre = Vector{FlexLevel}(undef, num_new_levels)
            for i in 1:num_new_levels
                u_bound = l_bound * 2.0
                pre[i] = FlexLevel((l_bound, u_bound), 0.0, 0.0, 0, Vector{Int64}())
                l_bound /= 2.0
            end
            prepend!(sampler.levels, pre)
        elseif extend_down
            num_new_levels = log_dist(l_bound, sampler.levels[end].bounds[1])
            post = Vector{FlexLevel}(undef, num_new_levels)
            for i in num_new_levels:-1:1
                u_bound = l_bound * 2.0
                post[i] = FlexLevel((l_bound, u_bound), 0.0, 0.0, 0, Vector{Int64}())
                l_bound = u_bound
            end
            append!(sampler.levels, post)
        else
            throw("sampler already contains FlexLevel of specified bounds")
        end
    end
    sampler.max_log2_upper_bound = floor_log2(sampler.levels[1].bounds[2])
end

"""
    trim_trailing_levels!(sampler)

Remove all empty `FlexLevel`s from the front and back of `sampler`.
"""
function trim_trailing_levels!(sampler::FlexleSampler)
    first = findfirst(level_is_populated, sampler.levels)
    if isnothing(first)
        sampler.levels = Vector{FlexLevel}()
        sampler.max_log2_upper_bound = nothing
    else
        last = findlast(level_is_populated, sampler.levels)
        sampler.levels = sampler.levels[first:last]
        sampler.max_log2_upper_bound = floor_log2(sampler.levels[1].bounds[2])
    end
end

"""
    print_flexle_sampler(sampler, name="")

Verbosely print `sampler` with optional label `name`.

Only used in testing; for end user printing, see [`Base.show`](@ref).
"""
function print_flexle_sampler(sampler::FlexleSampler; name::String="")
    @printf "\nFlexleSampler %s\n" name
    show(sampler.weights)
    @printf "\nSum %f" sampler.sum
    println()
    if isempty(sampler.levels)
        println("(sampler contains no levels)")
    else
        for l in sampler.levels
            @printf "\nLevel (%f, %f)\n" l.bounds[1] l.bounds[2]
            @printf "sum=%f\n" l.sum
            @printf "max=%f\n" l.max
            println(l.elements)
        end
    end
end

function Base.show(io::IO, sampler::FlexleSampler)
    l = length(sampler.weights)
    print(io, "FlexleSampler ($l weights)")
end

function print_weight(io::IO, sampler::FlexleSampler, i::Int64, l::Int64, pad::Int64)
    cap = i==l ? "" : "\n"
    print(io, lpad("$i", pad, " ") * ": $(sampler.weights[i])$cap")
end

function Base.show(io::IO, ::MIME"text/plain", sampler::FlexleSampler)
    l = length(sampler.weights)
    pad = 4 + ndigits(l)
    colon = iszero(l) ? "" : ":"
    print(io, "FlexleSampler ($l weights)$colon\n")
    n = numweights(sampler)
    if n > 20
        for i in 1:10
            print_weight(io, sampler, i, l, pad)
        end
        print(io, "    ⋮\n")
        for i in (n-8):n
            print_weight(io, sampler, i, l, pad)
        end
    else
        for i in eachindex(sampler.weights)
            print_weight(io, sampler, i, l, pad)
        end
    end
end
