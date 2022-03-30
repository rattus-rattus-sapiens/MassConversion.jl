using MassConversion
using Test

@testset verbose=true "MassConversion.jl" begin
    # Write your tests here.

    @testset "Data Structure Tests" begin
        traj1 = ones(10, 1000, 250)

        model = MassConversion.MCMOutput(traj1, 10.0, 0.1, [1.0, 1.0], [0.0], rand(Int64, 2, 2), [(0,0)], 100)
        @test size(model.trajectories_continuum) == (5, 1000, 250)
        @test size(model.trajectories_discrete)  == (5, 1000, 250)
        @test all(model.trajectories_continuum .== 1.0)
        @test all(model.trajectories_discrete .== 1.0)
    end

    @testset "Transition execution" begin
        state = [100, 100, 0, 100]
        MassConversion.exec_transition!(state, 2, 3)
        @test state == [100, 100, 0, 100]
        MassConversion.exec_transition!(state, 2, 1)
        @test state == [99, 100, 1, 100]
        MassConversion.exec_transition!(state, 2, 4)
        @test state == [99, 101, 1, 99]

        state = ones(100)
        MassConversion.exec_transition!(state, 50, 37)
        @test state[37] == 0
        @test state[87] == 2

        state[60] = 0.99999999
        MassConversion.exec_transition!(state, 50, 60)
        @test state[60] == 0
        @test state[10] == 2
    end

    @testset "Transition propensities" begin
        alpha = zeros(10)
        Kn = 5
        Rn = 0
        state = 50 .* ones(10)
        λ = ones(Float64, 5)
        θ = [(90, 110), (90, 90), (110, 110), (0, 0), (Inf, Inf)]
        enumθ = enumerate(θ)
        MassConversion.B!(alpha, state, λ, enumθ, Rn, Kn)
        @test alpha[1] == 0.0
        @test alpha[6] == 0.0
        @test alpha[2] == 50.0
        @test alpha[7] == 0.0
        @test alpha[3] == 0.0
        @test alpha[8] == 50.0
        @test alpha[4] == 50.0
        @test alpha[9] == 0.0
        @test alpha[5] == 0.0
        @test alpha[10] == 50.0
    end
end

#   @testset "Model setup tests" begin
#       model = Model()
#       model <= Species("A", 100)
#       model <= Species("B", 10)
#       model <= Species("EE", 9)
#
#       species = get_init(model)
#       @test species["A"] == 100
#       @test species["B"] == 10
#       @test species["EE"] == 9
#       @test_throws MassConversion.ModelSpecError Species("A+A", 0)
#       @test_throws MassConversion.ModelSpecError Species("A + 33", 4)
#       @test_throws MassConversion.ModelSpecError model <= Species("A", 10)
#    end
