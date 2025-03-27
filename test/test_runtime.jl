function testSampleRuntime!(w::Vector{Vector{Float64}}, w_sums::Vector{Float64}, sampler::FlexleSampler, r1::Int64, r2::Int64)
    r1 = 1:r1
    r2 = 1:r2
    
    println("CDF sampling:")    
    #randChoose(rand(), weights, w_tot)
    # @benchmark timeRandChoose($r1, $r2, $w, $w_tot)
    display(@benchmark timeRandChoose($r1, $r2, $w, $w_sums))
    # display(@benchmark randChoose(0.5, $weights, $w_tot, regenerate_rand=false))

    println("Flexle sampling:")
    # # sample(sampler)
    # @benchmark timeFlexleSample($r1, $r2, $sampler)
    display(@benchmark timeFlexleSample($r1, $r2, $sampler))
    # display(@benchmark(sample($sampler)))
end

function timeRandChoose(r1::UnitRange{Int64}, r2::UnitRange{Int64}, weights::Vector{Vector{Float64}}, w_sums::Vector{Float64})
    for i in r1
        rand_n = rand()
        for j in r2
            k = i*(j-1) + j
            y, rand_n = randChoose(rand_n, weights[k], w_sums[k], regenerate_rand=true)
        end
    end
end

function timeFlexleSample(r1::UnitRange{Int64}, r2::UnitRange{Int64}, sampler::FlexleSampler)
    for i in r1
        for j in r2
            y = sample(sampler)
        end
    end
end

"""
    testRuntime01(h=10000)

Compare the runtime of performing 1000 samples using CDF (i.e. `randChoose`, inverse transform) versus Flexle sampling on a
`Vector` of `h` weights, all either `0.0` or `1.0`.
"""
function testRuntime01(; h::Int64=10000, seed=0)
    Random.seed!(seed)
    w01::Vector{Float64} = zeros(h)
    for i in eachindex(w01)
        if rand() > 0.5
            w01[i] = 1
        end
    end

    r1, r2 = 200, 5
    w = Vector{Vector{Float64}}()
    w_sums = Vector{Float64}()
    for j in 1:(r1*r2)
        push!(w, zeros(h))
        for i in 1:h
            if rand() > 0.5
                w[j][i] = 1
            end
        end
        push!(w_sums, sum(w[j]))
    end


    testSampleRuntime!(w, w_sums, FlexleSampler(w01), r1, r2)
end

"""
    testRuntimeUniform(h=10000)

Compare the runtime of performing 1000 samples using CDF (i.e. `randChoose`, inverse transform) versus Flexle sampling on a
`Vector` of `h` weights, all randomly selected from a uniform distribution between `0.0` and `1.0`.
"""
function testRuntimeUniform(; h::Int64=10000, seed=0)
    Random.seed!(seed)
    wu::Vector{Float64} = rand(h)

    r1, r2 = 200, 5
    w = Vector{Vector{Float64}}()
    w_sums = Vector{Float64}()
    for j in 1:(r1*r2)
        push!(w, rand(h))
        push!(w_sums, sum(w[j]))
    end

    testSampleRuntime!(w, w_sums, FlexleSampler(wu), r1, r2)
end

function testStandardUpdateWeight!(weights::Vector{Float64}, updates::Tuple{Vector{Int64}, Vector{Float64}}, sum::Float64)
    i = updates[1]
    v = updates[2]
    for idx in eachindex(i)
        sum -= weights[i[idx]]
        weights[i[idx]] = v[idx]
        sum += v[idx]
    end
    return sum
end

function testFlexleUpdateWeight!(sampler::FlexleSampler, updates::Tuple{Vector{Int64}, Vector{Float64}})
    i = updates[1]
    v = updates[2]
    for idx in eachindex(i)
        sampler[i[idx]] = v[idx]
    end
    # verify(sampler)
end

"""
    testUpdateWeight(h=10000, n=1000, seed=0)

Compare the runtime of performing `n` random updates to a collection of `h` weights when storing said weights as
a simple `Vector` versus a `FlexleSampler`.
"""
function testUpdateWeight(; h::Int64=10000, n::Int64=1000, seed=0)
    Random.seed!(seed)
    w = rand(h)
    s = FlexleSampler(w)
    c = sum(w)
    updates = (rand(1:h, n), rand(n))

    println("Vector update weight:")    
    display(@benchmark testStandardUpdateWeight!($w, $updates, $c))

    println("Flexle update weight:")
    display(@benchmark testFlexleUpdateWeight!($s, $updates))
end

function testStandardAddWeight!(weights::Vector{Float64}, new_weights::Vector{Float64}, sum::Float64)
    for i in eachindex(new_weights)
        push!(weights, new_weights[i])
        sum += new_weights[i]
    end
    return sum
end

function testFlexleAddWeight!(sampler::FlexleSampler, new_weights::Vector{Float64})
    for i in eachindex(new_weights)
        push!(sampler, new_weights[i])
    end
    # verify(sampler)
end

"""
    testAddWeight(h=10000, n=1000, seed=0)

Compare the runtime of adding `n` random values to a starting collection of `h` weights when storing said weights as
a simple `Vector` versus a `FlexleSampler`.
"""
function testAddWeight(; h::Int64=10000, n::Int64=1000, seed=0)
    Random.seed!(seed)
    w = rand(h)
    s = FlexleSampler(w)
    c = sum(w)
    new_weights = rand(n)

    println("Vector add weight:")    
    display(@benchmark testStandardAddWeight!($w, $new_weights, $c))

    println("Flexle add weight:")
    display(@benchmark testFlexleAddWeight!($s, $new_weights))
end

function testStandardRemoveWeight!(weights::Vector{Float64}, indices::Vector{Int64}, sum::Float64)
    for i in eachindex(indices)
        sum -= weights[indices[i]]
        deleteat!(weights, indices[i])
    end
    return sum
end

function testFlexleRemoveWeight!(sampler::FlexleSampler, indices::Vector{Int64})
    for i in eachindex(indices)
        deleteat!(sampler, indices[i])
    end
    # verify(sampler)
end

"""
    testRemoveWeight(h=10000, n=1000, seed=0)

Compare the runtime of removing `n` random elements from a collection of `h` weights when storing said weights as
a simple `Vector` versus a `FlexleSampler`.
"""
function testRemoveWeight(; h::Int64=10000, n::Int64=1000, seed=0)
    Random.seed!(seed)
    w1 = rand(h)
    w2 = copy(w1)
    s = FlexleSampler(w2)
    c = sum(w1)
    indices = Vector{Int64}()
    for i in 1:n
        push!(indices, rand(1:h-i+1))
    end

    println("Vector remove weight:")    
    @time testStandardRemoveWeight!(w1, indices, c)

    println("Flexle remove weight:")
    @time testFlexleRemoveWeight!(s, indices)
end
