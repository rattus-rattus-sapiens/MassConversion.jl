using MassConversion
using Test

@testset "MassConversion.jl" begin
    # Write your tests here.
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
end
