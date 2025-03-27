using Flexle
using Test

include("flexletest.jl")


function runtests()
    @testset verbose=true "FlexleSampler user functions—creation/maintenance" begin
        @testset "FlexleSampler (initialize sampler from weights)" begin
            s = FlexleSampler(rand(100))
            @test verify(s, verbose=false) == 0

            s = FlexleSampler(zeros(100))
            @test verify(s, verbose=false) == 0

            s = FlexleSampler(zeros(100) .+ [((n >= 0.5) || (n < 0.25) ? n : 0.0) for n in rand(100)])
            @test verify(s, verbose=false) == 0

            s = FlexleSampler(ones(100) .* 1.5)
            @test verify(s, verbose=false) == 0

            s = FlexleSampler([10000.0 - 3.5*i for i in 1:1000])
            @test verify(s, verbose=false) == 0

            s = FlexleSampler(Vector{Float64}())
            @test verify(s, verbose=false) == 0
        end

        @testset "getindex (get existing weight)" begin
            s = FlexleSampler([100.0 - 2*i for i in 1:30])

            @test s[1] == 98.0
            @test s[2] == 96.0
            @test s[30] == 40.0

        end

        @testset "setindex! (update existing weights)" begin
            s = FlexleSampler([(i==27 ? 0.0 : 1.5 * i) for i in 1:100])

            # move between two levels - element in middle of indices vector
            s[24] = 24.0
            @test verify(s, verbose=false) == 0

            # move between two levels - element at end of indices vector
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
            s = FlexleSampler([rand()*10 for i in 1:1000])

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
            s = FlexleSampler([2.5, 6.0, 70.0, 0.001, 0.0, 4.2, 1.1])

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
        end
    end

end