module Flexle

include("sampler.jl")
include("interface.jl")

export FlexleSampler, update!, push!, deleteat!, getindex, setindex!, sample

end
