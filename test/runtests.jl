using BenchmarkTools
using Flexle
using HypothesisTests
using Printf
using ProfileView
using Test
using Random

include("test_correctness.jl")
include("test_runtime.jl")


@testset verbose=true "FlexleSampler initialization/maintenance" begin
    @testset "FlexleSampler (initialize sampler from weights)" begin
        s = FlexleSampler{Float64}(rand(100))
        @test verify(s, verbose=false) == 0

        s = FlexleSampler{Float64}(zeros(100))
        @test verify(s, verbose=false) == 0

        s = FlexleSampler{Float64}(zeros(100) .+ [((n >= 0.5) || (n < 0.25) ? n : 0.0) for n in rand(100)])
        @test verify(s, verbose=false) == 0

        s = FlexleSampler{Float64}(ones(100) .* 1.5)
        @test verify(s, verbose=false) == 0

        s = FlexleSampler{Float64}([10000.0 - 3.5*i for i in 1:1000])
        @test verify(s, verbose=false) == 0

        s = FlexleSampler{Float64}(Vector{Float64}())
        @test verify(s, verbose=false) == 0

        s = FlexleSampler{Int64}([5, 7, 1, 0, 0, 9, 0, 2, 10, 16, 400, 1, 0, 90])
        @test verify(s, verbose=false) == 0
    end

    @testset "getindex (get existing weight)" begin
        s = FlexleSampler{Float64}([100.0 - 2*i for i in 1:30])

        @test s[1] == 98.0
        @test s[2] == 96.0
        @test s[30] == 40.0

        s = FlexleSampler{Int64}([3*i for i in 1:15])
        @test s[3] == 9
    end

    @testset "setindex! (update existing weights)" begin
        s = FlexleSampler{Float64}([(i==27 ? 0.0 : 1.5 * i) for i in 1:100])

        # move between two levels - element in middle of elements vector
        s[24] = 24.0
        @test verify(s, verbose=false) == 0

        # move between two levels - element at end of elements vector
        s[85] = 60.0
        @test verify(s, verbose=false) == 0

        # within same level
        s[96] = 140.0
        @test verify(s, verbose=false) == 0

        # move between two levels - "from" level is now empty
        s[1] = 3.2
        @test verify(s, verbose=false) == 0

        # move between two levels - new level on top
        s[3] = 1000.0
        @test verify(s, verbose=false) == 0

        # move between two levels - new level on bottom
        s[40] = 0.001
        @test verify(s, verbose=false) == 0

        # make non-zero element zero
        s[50] = 0.0
        @test verify(s, verbose=false) == 0

        # make zero element non-zero, place in existing level
        s[27] = 500.1
        @test verify(s, verbose=false) == 0

        s[50] = 100000.0
        @test verify(s, verbose=false) == 0
    end
    
    @testset "push! (add new weight)" begin
        s = FlexleSampler{Float64}([rand()*10 for i in 1:1000])

        l = push!(s, 110.0)
        @test (verify(s, verbose=false) == 0) && (s[length(s.weights)] == 110.0)

        push!(s, 80.0)
        @test (verify(s, verbose=false) == 0) && (s[length(s.weights)] == 80.0)

        push!(s, 0.0)
        @test (verify(s, verbose=false) == 0) && (s[length(s.weights)] == 0.0)

        push!(s, 0.126)
        @test (verify(s, verbose=false) == 0) && (s[length(s.weights)] == 0.126)

        push!(s, 1.0)
        @test (verify(s, verbose=false) == 0) && (s[length(s.weights)] == 1.0)
    end

    @testset "deleteat! (remove existing weight)" begin
        s = FlexleSampler{Float64}([2.5, 6.0, 70.0, 0.001, 0.0, 4.2, 1.1])

        deleteat!(s, 2)     # 6.0, remove from level w/o emptying
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 6)

        deleteat!(s, 5)     # 4.2, remove from middle level and empty
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 5)

        deleteat!(s, 2)     # 70.0, empty top level
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 4)

        deleteat!(s, 3)     # 0.0, remove zero weight
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 3)

        deleteat!(s, 2)     # 0.001, empty bottom level
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 2)

        deleteat!(s, 1)     # 2.5, remove first element
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 1)

        deleteat!(s, 1)     # 1.1, empty sampler
        @test (verify(s, verbose=false) == 0) && (length(s.weights) == 0)

        s = FlexleSampler{Float64}([0.0, 2.0, 0.0, 1.0, 1.0, 0.0, 2.0, 2.0, 1.0, 0.0, 0.0, 1.0, 1.0])
        for i in 1:numweights(s)
            deleteat!(s, 1)
            @test (verify(s, verbose=false) == 0)
        end
    end
end

@testset "Flexle stats" begin
    s = FlexleSampler{Float64}([2.0*i for i in 1:5])

    @test s[2] == 4.0
    @test s[5] == 10.0
    @test_throws BoundsError s[6]

    @test getweights(s) == [2.0, 4.0, 6.0, 8.0, 10.0]
    @test numweights(s) == 5

    push!(s, 11.0)
    push!(s, 12.0)
    deleteat!(s, 1)
    @test getweights(s) == [4.0, 6.0, 8.0, 10.0, 11.0, 12.0]
    @test numweights(s) == 6
    @test_throws BoundsError s[0]
end

@testset "Flexle sampling distribution" begin
    alpha = 0.01

    # test for empirical distribution _not_ significantly different from expected given weights.
    # each p-value tests will fail ~alpha% of the time by chance.
    # correspondingly, the chance of _no_ tests failing is (1-alpha)^num_tests (~96% for 4 @ alpha=0.01)

    s = FlexleSampler{Float64}([rand() for _ in 1:20])
    @test pvalue(flexleChiSquared(s)) > alpha

    s = FlexleSampler{Float64}([20.0-i for i in 1:20])
    @test pvalue(flexleChiSquared(s)) > alpha

    push!(s, 6.7)
    push!(s, 0.0)
    deleteat!(s, 2)
    s[3] = 10.1
    @test pvalue(flexleChiSquared(s)) > alpha

    s = FlexleSampler{Float64}(vcat(zeros(10), [1.0], zeros(5), [2.0]))
    @test pvalue(flexleChiSquared(s)) > alpha

    s = FlexleSampler{Float64}(zeros(1000))
    @test_throws DomainError flexleChiSquared(s)    # can't sample from all-0 distribution


end
