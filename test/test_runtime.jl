using StatsBase
using Plots

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



"""
Benchmark tests for README.md.
"""
const BLUE::Int64 = 1
const GREEN::Int64 = 3
const PURPLE::Int64 = 4


function compare_sampling()
    sizes = [10, 100, 1000, 10000, 100000, 1000000]
    weight_vectors = [rand(n) for n in sizes]

    weights  = [Weights(w) for w in weight_vectors]
    samplers = [FlexleSampler(w) for w in weight_vectors]
    
    weights_results::Vector{Float64} = [mean(@benchmark(StatsBase.sample($w))).time for w in weights]
    flexle_results::Vector{Float64} = [mean(@benchmark(Flexle.sample($s))).time for s in samplers]

    return sizes, weights_results, flexle_results
end

function plot_compare_sampling(; path="docs/assets/", extension=".png", write=true)
    sizes, weights_results, flexle_results = compare_sampling()
    num_groups = length(sizes)
    num_categories = 2

    ctg = repeat(["StatsBase default sample", "Flexle sample"], inner = num_groups)
    nam = repeat([string(s) for s in sizes], outer = num_categories)

    
    p = groupedbar(nam, hcat(weights_results, flexle_results), group = ctg, xlabel = "#weights", ylabel = "Mean sample time (ns)",
        # title = "Scores by group and category", bar_width = 0.67,
        lw = 0, yaxis=:log, ylims=(1e0, 1e+6), framestyle = :box, legend=:topleft, color=repeat([BLUE, GREEN], inner=num_groups))

    if write
        savefig(p, path * "01_compare_sampling" * extension)
    end

    display(p)
end

function alias_runtime_test(weights, n_samples)
    indices = eachindex(weights)
    result = Vector{Int64}(undef, n_samples)
    StatsBase.alias_sample!(indices, Weights(weights), result)
    return result
end

function flexle_runtime_test(weights, n_samples)
    sampler = FlexleSampler(weights)
    result = [Flexle.sample(sampler) for _ in 1:n_samples]
    return result
end

function compare_sampling_alias(num_samples::Vector{Int64}, sizes::Vector{Int64})
    weight_vectors = [rand(n) for n in sizes]

    alias_results  = Matrix{Float64}(undef, length(num_samples), length(sizes))
    flexle_results = Matrix{Float64}(undef, length(num_samples), length(sizes))

    for j in eachindex(sizes)
        for i in eachindex(num_samples)
            alias_results[i,j] = mean(@benchmark(alias_runtime_test($(weight_vectors[j]), $(num_samples[i])))).time
            flexle_results[i,j] = mean(@benchmark(flexle_runtime_test($(weight_vectors[j]), $(num_samples[i])))).time
            println("ran: 10^", i, " samples, 10^", j, " weights")
        end
    end

    return alias_results, flexle_results
end

# Warning: this function takes a while to run (on the order of ~ten minutes on my Intel i7, 8GB memory laptop)
function plot_compare_sampling_alias(; path="docs/assets/", extension=".png", write=true)
    sample_size_range = 1:6
    weightvector_size_range = 1:6

    num_samples = [10^k for k in sample_size_range]
    sizes       = [10^l for l in weightvector_size_range]

    alias_results, flexle_results = compare_sampling_alias(num_samples, sizes)
    runtime_ratio_results = log2.(flexle_results ./ alias_results)

    lim = ceil(maximum(abs.(runtime_ratio_results)))    # center 0 in heatmap by making upper and lower limits equal
    p = heatmap([string(l) for l in weightvector_size_range],
                [string(k) for k in sample_size_range],
                runtime_ratio_results, xlabel = "log_10(#weights)", ylabel = "log_10(#samples)", clim=(-lim, lim), seriescolor=:bluesreds)

    if write
        savefig(p, path * "02_compare_sampling_alias" * extension)
    end

    display(p)
end

function first_decimal_place(n::Float64)
    return Int64(round(10.0 * (trunc(n, digits=1) - floor(n))))
end


function do_ops_statsbase(weights::Vector{Float64}, r::Int64)
    w = Weights(weights)
    for _ in 1:r
        i = StatsBase.sample(w)
        v = w[i]
        f = first_decimal_place(v)

        if f == 0
            w = Weights(deleteat!(w.values, i))
        elseif f == 5
            w = Weights(push!(w.values, v/7.0))
        else
            w[i] = 2.0*v
        end
    end
    return w
end

function do_ops_flexle(weights::Vector{Float64}, r::Int64)
    sampler = FlexleSampler(weights)
    for _ in 1:r
        i = Flexle.sample(sampler)
        v = sampler[i]
        f = first_decimal_place(v)

        if f == 0
            deleteat!(sampler, i)
        elseif f == 5
            push!(sampler, v/7.0)
        else
            sampler[i] = 2.0*v
        end
    end
    return sampler
end

function compare_ops(num_iterations::Vector{Int64}, sizes::Vector{Int64})
    weight_vectors = [rand(n) for n in sizes]

    statsbase_results  = Matrix{Float64}(undef, length(num_iterations), length(sizes))
    flexle_results = Matrix{Float64}(undef, length(num_iterations), length(sizes))

    for j in eachindex(sizes)
        for i in eachindex(num_iterations)
            statsbase_results[i,j] = mean(@benchmark(do_ops_statsbase($(weight_vectors[j]), $(num_iterations[i])))).time
            flexle_results[i,j] = mean(@benchmark(do_ops_flexle($(weight_vectors[j]), $(num_iterations[i])))).time
            println("ran: 10^", i, " iterations, 10^", j, " weights")
        end
    end

    return statsbase_results, flexle_results
end

function plot_compare_ops(; path="docs/assets/", extension=".png", write=true)
    Random.seed!(1)
    num_iteration_range = 1:5
    weightvector_size_range = 1:5

    num_iterations = [10^k for k in num_iteration_range]
    sizes          = [10^l for l in weightvector_size_range]

    statsbase_results, flexle_results = compare_ops(num_iterations, sizes)
    runtime_ratio_results = log2.(flexle_results ./ statsbase_results)

    lim = ceil(maximum(abs.(runtime_ratio_results)))    # center 0 in heatmap by making upper and lower limits equal
    p = heatmap([string(l) for l in weightvector_size_range],
                [string(k) for k in num_iteration_range],
                runtime_ratio_results, xlabel = "log_10(#weights)", ylabel = "log_10(#iterations)", clim=(-lim, lim), seriescolor=:bluesreds)

    if write
        savefig(p, path * "03_compare_ops" * extension)
    end

    display(p)
end
