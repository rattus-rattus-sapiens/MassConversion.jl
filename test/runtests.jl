using MassConversion
using Test
using StaticArrays

@testset verbose=true "MassConversion.jl" begin

    t_span       = 5.0:0.01:10
    D_init       = rand(Int64, 3)
    C_init       = rand(Float64, 3)
    λ_reac       = rand(5)
    λ_tran       = rand()
    R_mat        = rand(Int64, 6, 5)
    θ            = [(0,0), (0,0), (0,0)]
    calc_prop(x) = x
    calc_dxdt(x) = x

    model  = MassConversion.MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mat, θ, calc_prop, calc_dxdt)
    state  = MassConversion.MCMstate(model)
    output = MassConversion.MCMraw(model)

    @testset verbose=true "MCMmodel Tests" begin

        @testset "scalar params" begin
            @test model.t_start == 5.0
            @test model.t_step  == 0.01
            @test model.t_final == 10.0
            @test model.t_len   == 501
            @test model.n_spec  == 3
            @test model.n_reac  == 5
        end

        @testset "arrays" begin
            @test model.D_init == D_init
            @test model.C_init == C_init
            @test model.R_mat  == R_mat
            @test model.λ_reac == λ_reac
            @test model.λ_tran == repeat([λ_tran], 6)
            @test model.θ      == θ

            @test model.D_init isa SVector{3, Int64}
            @test model.C_init isa SVector{3, Float64}
            @test model.R_mat  isa Matrix{Int64}
            @test model.λ_reac isa SVector{5, Float64}
            @test model.λ_tran isa SVector{6, Float64}
            @test model.θ       isa Vector{Tuple{Float64, Float64}}
        end
    end

    @testset verbose=true "MCMraw Tests" begin
        @testset "scalar params" begin
            @test output.n == 0
        end
    end
end