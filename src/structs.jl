# Flexle structs

const WeightNumber = Union{Float64, Int64}  # TODO: any Int or Float

mutable struct FlexLevel{T <: WeightNumber}
    bounds::Tuple{T,T}  # lower, upper
    sum::T
    max::T
    num_max::Int64
    elements::Vector{Int64}
end

mutable struct FlexleSampler{T <: WeightNumber}
    levels::Vector{FlexLevel{T}}
    weights::Vector{T}
    sum::T
    element_positions::Vector{Int64}
    max_log2_upper_bound::Union{Int64, Nothing}     # nothing when levels is empty
end
