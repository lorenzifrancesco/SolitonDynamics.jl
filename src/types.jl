abstract type Params end
abstract type TransformLibrary end
abstract type Space end
# abstract parameter type: can also be a function of simulation time

abstract type PotentialType end
abstract type Gaussian <: PotentialType end
abstract type Elliptical <: PotentialType end

struct EquationType
    D::Int64
    improved::Bool
end

const GPE_1D = EquationType(1, false)
const NPSE = EquationType(1, true)
const GPE_3D = EquationType(3, false)

abstract type Simulation{D} end

abstract type Solver end
abstract type SplitStep <: Solver end
abstract type CrankNicholson <: Solver end
abstract type UserParams end
abstract type Method end

struct XSpace{D} <: Space
    psiX::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

struct KSpace{D} <: Space
    psiK::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

function crandn_array(D, N, t)
    if t == CuArray
        return 
    else
    return rand(N...)
end

function Params() <: UserParams
    return 0.0
end

@with_kw mutable struct Transforms{D,N} <: TransformLibrary
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D,UnitRange{Int64}},Float64} = 0.1* plan_fft(crandn_array(D, N, t))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D,UnitRange{Int64}},Float64} = 0.1*plan_fft!(crandn_array(D, N, t))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D,UnitRange{Int64}},Float64} = 0.1* plan_ifft(crandn_array(D, N, t))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D,UnitRange{Int64}},Float64} = 0.1*plan_ifft!(crandn_array(D, N, t))
    #psi::ArrayPartition = crandnpartition(D,N)
end

@with_kw mutable struct Sim{D} <: Simulation{D}
    # === solver and algorithms
    equation::EquationType
    solver::Solver = SplitStep
    graphics::Bool = false
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm = Tsit5() # default solver
    reltol::Float64 = 1e-6 # default tolerance; may need to use 1e-7 for corner cases
    flags::UInt32 = FFTW.MEASURE # choose a plan. PATIENT, NO_TIMELIMIT, EXHAUSTIVE
    # === dimensions and physics
    L::NTuple{D,Float64} # length scales
    N::NTuple{D,Int64}  # grid points in each dimensions
    g = 0.1
    γ = 0.0; @assert γ >= 0.0 # three body losses param
    ti = 0.0    # initial time
    tf = 2    # final time
    Nt::Int64 = 200     # number of saves over (ti,tf)
    params::UserParams = Params() # optional user parameterss
    V0::AbstractArray{Float64,D}
    t::LinRange{Float64} = LinRange(ti,tf,Nt) # time of saves
    psi_0::AbstractArray{Complex{Float64},D} # initial condition

    # === saving
    nfiles::Bool = false
    path::String = nfiles ? joinpath(@__DIR__,"data") : @__DIR__
    filename::String = "save"
    # === arrays, transforms, spectral operators
    X::NTuple{D,AbstractArray{Float64,1}} = xvecs(L,N)
    K::NTuple{D,AbstractArray{Float64,1}} = kvecs(L,N)
    T::TransformLibrary = makeT(X,K,flags=flags)
end

InitSim(L,N,par) = Sim{length(L)}(L=L,N=N,params=par)
InitSim(L,N) = Sim{length(L)}(L=L,N=N,params=Params())