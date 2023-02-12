abstract type Simulation{D} end
abstract type Params end
abstract type Transforms end
abstract type Field end
# abstract parameter type: can also be a function of simulation time

struct Xfield{D} <: Field
    psiX::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

struct Kfield{D} <: Field
    psiK::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

@with_kw mutable struct Sim{D} <: Simulation{D}
    k 
end