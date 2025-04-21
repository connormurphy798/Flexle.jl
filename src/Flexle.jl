module Flexle

include("structs.jl")
include("math.jl")
include("maintenance.jl")
include("sampling.jl")
include("userinterface.jl")

export FlexleSampler, getindex, setindex!, getweights, numweights, push!, deleteat!, sample

end
