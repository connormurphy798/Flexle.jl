# Flexle structs

mutable struct FlexLevel
    bounds::Tuple{Float64,Float64}  # lower, upper
    sum::Float64
    max::Float64
    num_max::Int64
    elements::Vector{Int64}
end

mutable struct FlexleSampler
    levels::Vector{FlexLevel}
    weights::Vector{Float64}
    sum::Float64
    element_positions::Vector{Int64}
    max_log2_upper_bound::Union{Int64, Nothing}     # nothing when levels is empty
end
