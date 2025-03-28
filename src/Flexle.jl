module Flexle

include("sampler.jl")
include("interface.jl")

export FlexleSampler, getindex, setindex!, getweights, numweights, push!, deleteat!, sample

end
