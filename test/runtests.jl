using SolitonDynamics, Test

@testset "Transform" begin
    include("test_transforms.jl")
end
@testset "Ground state" begin
    include("test_ground_state.jl")
end
@testset "Dynamics" begin
    include("test_dynamics.jl")
end
