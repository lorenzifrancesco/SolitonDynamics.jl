abstract type TransformLibrary{A <: AbstractArray} end
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

struct Solver
    spectral::Bool;
end

const SplitStep = Solver(true)
const CrankNicholson = Solver(false)

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

@with_kw mutable struct Params <: UserParams
    κ = 0.0 # a placeholder
end

@with_kw mutable struct Transforms{D,N,A} <: TransformLibrary{A}
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D,UnitRange{Int64}},Float64} = 0.1* plan_fft(crandn_array(N, A))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D,UnitRange{Int64}},Float64} = 0.1*plan_fft!(crandn_array(N, A))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D,UnitRange{Int64}},Float64} = 0.1* plan_ifft(crandn_array(N, A))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D,UnitRange{Int64}},Float64} = 0.1*plan_ifft!(crandn_array(N, A))
    #psi::ArrayPartition = crandnpartition(D,N,A)
end

@with_kw mutable struct Sim{D, A <: AbstractArray}
    # === solver and algorithm
    equation::EquationType = GPE_1D
    solver::Solver = SplitStep
    graphics::Bool = false
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm = Tsit5() # default solver
    reltol::Float64 = 1e-6 # default tolerance; may need to use 1e-7 for corner cases
    flags::UInt32 = FFTW.MEASURE # choose a plan. PATIENT, NO_TIMELIMIT, EXHAUSTIVE
    iswitch::ComplexF64 = 1.0 # 1.0 for real time, -im for imaginary time
    # === dimensions and physics
    L::NTuple{D,Float64} # length scales
    N::NTuple{D,Int64}  # grid points in each dimensions
    g = 0.1
    γ = 0.0; @assert γ >= 0.0 # three body losses param
    
    ti = 0.0    # initial time
    tf = 1.0    # final time
    tspan = [ti, tf]
    Nt::Int64 = 5     # number of saves over (ti,tf)
    params::UserParams = Params() # optional user parameterss
    V0::A = zeros(N)
    t::LinRange{Float64} = LinRange(ti,tf,Nt) # time of saves
    psi_0::A = zeros(N) |> complex # initial condition
    dV = volume_element(L, N)
    Vol = prod(L)
    # === saving
    nfiles::Bool = false
    path::String = nfiles ? joinpath(@__DIR__,"data") : @__DIR__
    filename::String = "save"
    # === arrays, transforms, spectral operators
    X::NTuple{D,A} = xvecs(L,N)
    K::NTuple{D,A} = kvecs(L,N)
    T::TransformLibrary{A} = makeT(X,K,A,flags=flags)
    ksquared::A = 0.5*k2(K, A)
end

InitSim(L,N,A,par) = Sim{length(L), A}(L=L,N=N,params=par)
InitSim(L,N,A) = Sim{length(L), A}(L=L,N=N,params=Params())