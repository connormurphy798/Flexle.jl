using BenchmarkTools
using Printf
using Random
using ProfileView

"""
    verify(sampler, name="")

Performs several checks on `sampler` to ensure that it is internally consistent, returning the total number
of inconsistencies detected.

`name` allows the caller to give an optional label to `sampler` for pretty printing purposes.

Checks are as follows:
- index placement in `sampler.levels`
    - every index of non-zero weight in `sample.weights` appears in some level
    - every index appearing in a level is an index that:
        - exists in `sampler.weights`, and if so,
        - belongs in this level according to its weight
- construction of `index_positions` at each level in `sampler.levels`
    - each index appears in the appropriate level's `index_positions` corresponding to its position in `level.indices`
    - no index appears in an `index_positions` corresponding to a position at which it does not exist
- `sampler` stats
    - each level's `sum` is (approximately) equal to the sum of the weights of the indexes it holds
    - each level's `max` is equal to the maximum of the weights of the indexes it holds
    - `sampler.sum` is (approximately) equal to the sum of `sampler.weights`
    - `sampler.max_log2_upper_bound` is equal to the log2 of the highest level's upper bound
"""
function verify(sampler::FlexleSampler; name::String="")
    @printf "Verifying FlexleSampler %s...\n" name
    errors = 0

    # all weights in sampler.weights are represented in some level?
    for i in eachindex(sampler.weights)
        if !iszero(sampler.weights[i]) && !Flexle.inSampler(i, sampler)
            errors += 1
            @printf "Error %i: index %i (weight %f) not in any level\n" errors i sampler.weights[i]
        end
    end

    # all weights in levels are represented in sampler.weights, and if so, in the correct level?
    num_weights = length(sampler.weights)
    for level in sampler.levels
        for i in level.indices
            if !(1 <= i <= num_weights)
                errors += 1
                @printf "Error %i: index %i present in level (%f, %f) but not in weights\n" errors i level.bounds[1] level.bounds[2]
            elseif !(level.bounds[1] <= sampler.weights[i] < level.bounds[2])
                errors += 1
                @printf "Error %i: index %i (weight %f) present in level (%f, %f)\n" errors i sampler.weights[i] level.bounds[1] level.bounds[2]
            end
        end
    end

    # all indices in levels recorded correctedly in index_positions?
    d = sampler.index_positions
    for level in sampler.levels
        for pos in eachindex(level.indices)
            idx = level.indices[pos]
            if d[idx] != pos
                errors += 1
                @printf "Error %i: element %i (level (%f, %f), weight %f) not recorded as position %i in index_positions\n" errors idx level.bounds[1] level.bounds[2] sampler.weights[idx] pos
            end
        end
    end
    for idx in eachindex(d)
        pos = d[idx]
        level = Flexle.getLevel(sampler.weights[idx], sampler)
        if !iszero(pos) && ((pos > length(level.indices)) || !(level.indices[pos] == idx))
            errors += 1
            @printf "Error %i: element %i (level (%f, %f), weight %f) incorrectly recorded as position %i in index_positions\n" errors idx level.bounds[1] level.bounds[2] sampler.weights[idx] pos
        end
    end

    # all recorded stats correct?
    overall_sum = 0.0
    for level in sampler.levels
        empty = isempty(level.indices)
        expected_sum = empty ? 0.0 : sum(sampler.weights[i] for i in level.indices)
        if !Flexle.approxeq(expected_sum, level.sum)
            errors += 1
            @printf "Error %i: level (%f, %f) sum incorrect (expected %f, got %f)\n" errors level.bounds[1] level.bounds[2] expected_sum level.sum
        end
        overall_sum += expected_sum

        expected_max = empty ? 0.0 : maximum(sampler.weights[i] for i in level.indices)
        if expected_max != level.max
            errors += 1
            @printf "Error %i: level (%f, %f) max incorrect (expected %f, got %f)\n" errors level.bounds[1] level.bounds[2] expected_max level.max
        end
    end
    if !Flexle.approxeq(overall_sum, sampler.sum)     # correct for probable floating point error
        errors += 1
        @printf "Error %i: overall sampler sum incorrect (expected %f, got %f)" errors overall_sum sampler.sum
    end
    if !((isempty(sampler.levels) && isnothing(sampler.max_log2_upper_bound)) || (!isempty(sampler.levels) && (sampler.max_log2_upper_bound == Flexle.floorLog2(sampler.levels[1].bounds[2]))))
        errors += 1
        @printf "Error %i: sampler max_log2_upper_bound incorrect (expected %i, got %i)" errors Flexle.floorLog2(sampler.levels[1].bounds[2]) sampler.max_log2_upper_bound
    end

    @printf"Final error count: %i\n" errors
    return errors
end

function randChoose(rand_n::Float64, rates::AbstractVector{Float64}, rates_sum::Float64; regenerate_rand::Bool=false)
    #TODO: rates is not type-defined to accommodate @views, maybe this is bad
    if rates_sum > 0.0
        norm_rand_n = rand_n * rates_sum
        cum_sum = 0.0
        idx = 0
        chosen_rate = 0.0
        for r in rates
            idx += 1
            cum_sum += r
            chosen_rate = r
            if cum_sum > norm_rand_n
                break
            end
        end

        if regenerate_rand
            return idx, (norm_rand_n - cum_sum + chosen_rate) / chosen_rate
        else
            return idx, 0.0
        end
    else
        return NaN, rand_n
    end
end

"""
    testFlexleSampler()

Test the basic internal maintenance of the `FlexleSampler` struct.

Assesses:
- initialization
- updating of weights
- extending of levels
"""
function testFlexleSampler(seed=0)
    Random.seed!(seed)

    # --- INITIALIZATION ---
    println("\nTEST: FlexleSampler initialization from weights vectors\n")
    names = Vector{String}()
    weightss = Vector{Vector{Float64}}()
    samplers = Vector{FlexleSampler}()

    push!(names, "(random weights vector, len 100)")
    push!(weightss, rand(100))
    push!(samplers, FlexleSampler(weightss[end]))
    verify(samplers[end], name=names[end])

    push!(names, "(all zero weights vector, len 100)")
    push!(weightss, zeros(100))
    push!(samplers, FlexleSampler(weightss[end]))
    verify(samplers[end], name=names[end])

    push!(names, "(some zeros in weights vector, len 100)")
    push!(weightss, zeros(100) .+ [((n >= 0.5) || (n < 0.25) ? n : 0.0) for n in rand(100)])
    push!(samplers, FlexleSampler(weightss[end]))
    verify(samplers[end], name=names[end])


    # --- MOVE BETWEEN EXISTING LEVELS ---

    push!(names, "(move between levels, test 1)")
    push!(weightss, [1.5 * n for n in 1.0:100.0])
    push!(samplers, FlexleSampler(weightss[end]))
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(move test, pre-move)")

    # between two levels - element in middle of indices vector
    idx = 24
    w_new = 24.0 # now belongs in (16.0, 32.0) instead of (32.0, 64.0)
    samplers[end][idx] = w_new
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(move test, changed weights[24] from 36.0 to 24.0)")

    # between two levels - element at end of indices vector
    idx = 85
    w_new = 60.0 # now belongs in (32.0, 64.0) instead of (64.0, 128.0)
    samplers[end][idx] = w_new
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(move test, changed weights[85] from 127.5 to 60.0)")

    # within same level
    idx = 96
    w_new = 140.0 # still belongs in (128.0, 256.0)
    samplers[end][idx] = w_new
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(move test, changed weights[96] from 144.0 to 140.0)")

    # between two levels - "from" level is now empty
    idx = 1
    w_new = 3.2 # now belongs in (2.0, 4.0) instead of (1.0, 2.0)
    samplers[end][idx] = w_new
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(move test, changed weights[1] from 1.5 to 3.2)")


    # --- EXTEND LEVELS ---

    # upwards by 1
    Flexle.extendLevels!((256.0, 512.0), samplers[end])
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, single level upwards)")

    # upwards by several
    Flexle.extendLevels!((4096.0, 8192.0), samplers[end])
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, several levels upwards)")

    # downwards by one
    Flexle.extendLevels!((0.5, 1.0), samplers[end])
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, single level downwards)")

    # downwards by several
    Flexle.extendLevels!((0.015625, 0.03125), samplers[end])
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, several levels downards)")

    # ERROR: new bound matches existing upper bound
    try
        Flexle.extendLevels!((4096.0, 8192.0), samplers[end])
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, attempt to add existing upper bound)")

    #  ERROR: new bound matches existing lower bound
    try
        Flexle.extendLevels!((0.015625, 0.03125), samplers[end])
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, attempt to add existing lower bound)")

    # ERROR: new bound matches occupied intermediate level
    try
        Flexle.extendLevels!((2.0, 4.0), samplers[end])
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, attempt to add existing occupied intermediate level)")

    # ERROR: new bound matches empty intermediate level
    try
        Flexle.extendLevels!((1.0, 2.0), samplers[end])
    catch e
        println("correctly identified error: ", e)
    end
    # printFlexleSampler(samplers[end])
    verify(samplers[end], name="(extend test, attempt to add existing empty intermediate level)")
end

"""
    levelSamplingExpectedValue(sampler, n=1)

Produce a `Dict` mapping (bounds of) levels in `sampler` to their expected number of selections given `n` samples.

For default value `n=1`, each level is simply assigned the probability of its selection.
"""
function levelSamplingExpectedValue(sampler::FlexleSampler; n::Int64=1)
    norm = Float64(n) / sampler.sum
    d = Dict{Tuple{Float64,Float64},Float64}()
    for level in sampler.levels
        d[level.bounds] = norm * level.sum
    end
    return d
end

"""
    testLevelSampling(sampler, n)

Randomly select a level from `sampler` `n` times and return a `Dict` mapping levels to empirical number of draws.
"""
function testLevelSampling(sampler::FlexleSampler, n::Int64)
    d = Dict{Tuple{Float64,Float64},Int64}()
    for level in sampler.levels
        d[level.bounds] = 0
    end

    for _ in 1:n
        bounds = Flexle.cdfSample(sampler)[1].bounds
        d[bounds] += 1
    end

    return d
end

"""
    indexSamplingExpectedValue(level, weights, n=1)

Produce a `Dict` mapping indexes in `level` to their expected number of selections (according to `weights`) given `n` samples.

For default value `n=1`, each index is simply assigned the probability of its selection. 
"""
function indexSamplingExpectedValue(level::Flexle.FlexLevel, weights::Vector{Float64}; n::Int64=1)
    norm = Float64(n) / level.sum
    d = Dict{Int64,Float64}()
    for i in level.indices
        d[i] = norm * weights[i]
    end
    return d
end

"""
    testIndexSampling(level, weights, n)

Randomly select an index from `level` `n` times according to `weights` and return a `Dict` mapping indexes to empirical number of draws.
"""
function testIndexSampling(level::Flexle.FlexLevel, weights::Vector{Float64}, n::Int64)
    d = Dict{Int64,Int64}()
    for i in level.indices
        d[i] = 0
    end

    for _ in 1:n
        i = Flexle.rejectionSample(rand(), level, weights)
        d[i] += 1
    end

    return d
end

"""
    flexleSamplingExpectedValue(sampler, n=1)

Produce a `Dict` mapping indexes in `sampler` to their expected number of selections given `n` samples.

For default value `n=1`, each index is simply assigned the probability of its selection. 
"""
function flexleSamplingExpectedValue(sampler::FlexleSampler; n::Int64=1)
    norm = Float64(n) / sampler.sum # (l.sum / sampler.sum) * (Float64(n) / l.sum) -- l terms cancel out
    d = Dict{Int64,Float64}()
    for i in eachindex(sampler.weights)
        d[i] = norm * sampler.weights[i]
    end

    # for l in sampler.levels
    #     for i in l.indices
    #         d[i] = norm * sampler.weights[i]
    #     end
    # end
    return d
end

"""
    testFlexleSampling(sampler, n)

Randomly select an index from `sampler` `n` times and return a `Dict` mapping indexes to empirical number of draws.
"""
function testFlexleSampling(sampler::FlexleSampler, n::Int64)
    d = Dict{Int64,Int64}()
    for i in eachindex(sampler.weights)
        d[i] = 0
    end

    for _ in 1:n
        i = sample(sampler)
        d[i] += 1
    end

    return d
end

"""
    testUserFunctions()

Test the internal maintenance of the `FlexleSampler` struct in response to user commands.

Assesses:
- initialization
- updating of weights
- addition of new weights
- removal of existing weights
"""
function testUserFunctions()
    # test create
    w0 = zeros(1000)
    s0 = FlexleSampler(w0)
    verify(s0)

    w1 = w0 .+ [((i < 0.25 || 0.50 < i) ? i : 0.0) for i in rand(1000)]
    s1 = FlexleSampler(w1)
    verify(s1)

    w2 = zeros(100) .+ [(i > 0.5 ? i : 0.0) for i in 0.0:99.0]
    s2 = FlexleSampler(w2)
    verify(s2)

    # test update
    for i in 2:10:92
        s2[i] = 100*rand()
    end
    verify(s2)
    
    s2[87] = 0.0     # make non-zero element 0.0
    verify(s2)

    s2[1] = 3000.0   # make 0.0 element non-zero
    verify(s2)

    # test add
    push!(s2, 99.5)       # in existing level (occupied)
    push!(s2, 480.0)      # in existing level (empty)
    push!(s2, 10000.0)    # in non-existing level (above)
    push!(s2, 0.001)      # in non-existing level (below)
    push!(s2, 0.0)        # zero
    verify(s2)

    # test remove
    deleteat!(s2, 70)      # in middle level, does not empty level
    deleteat!(s2, 101)     # in middle level, empties level
    deleteat!(s2, 101)     # empty top level
    deleteat!(s2, 101)     # empty bottom level
    deleteat!(s2, 86)      # zero
    deleteat!(s2, 1)       # first element
    verify(s2)
end

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
    flexleTestSuite()

Calls a series of Flexle test function to ensure working status before deployment.
"""
function flexleTestSuite()
    # functionality
    testFlexleSampler()
    testUserFunctions()

    # runtime
    testRuntime01()
    testRuntimeUniform()
    testUpdateWeight()
    testAddWeight()
    testRemoveWeight()
end
