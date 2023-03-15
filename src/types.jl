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
    number::Int64
    spectral::Bool;
end

const SplitStep = Solver(1, true)
const CrankNicholson = Solver(2, false)
const PredictorCorrector = Solver(3, false)
const BackwardEuler = Solver(4, false)


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

@with_kw mutable struct GPUTransforms{D,N,A} <: TransformLibrary{A}
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64}, CUDA.CUFFT.cCuFFTPlan{ComplexF64, -1, false, 3},Float64} = 0.1*CUDA.CUFFT.plan_fft(  CuArray(crandn_array(N, A)))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},CUDA.CUFFT.cCuFFTPlan{ComplexF64, -1, true , 3},Float64} = 0.1*CUDA.CUFFT.plan_fft!( CuArray(crandn_array(N, A)))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64}, CUDA.CUFFT.cCuFFTPlan{ComplexF64,  1, false, 3},Float64} = 0.1*CUDA.CUFFT.plan_ifft( CuArray(crandn_array(N, A)))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},CUDA.CUFFT.cCuFFTPlan{ComplexF64,  1, true , 3},Float64} = 0.1*CUDA.CUFFT.plan_ifft!(CuArray(crandn_array(N, A)))
    #psi::ArrayPartition = crandnpartition(D,N,A)
end

@with_kw mutable struct Sim{D, A <: AbstractArray}

    # === solver and algorithm
    equation::EquationType = GPE_1D
    solver::Solver = SplitStep
    graphics::Bool = false
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm = Tsit5() # default solver
    reltol::Float64 = 1e-3  # default tolerance; may need to use 1e-7 for corner cases
    abstol::Float64 = 1e-3
    maxiters::Int64 = 5000
    flags::UInt32 = FFTW.MEASURE # choose a plan. PATIENT, NO_TIMELIMIT, EXHAUSTIVE
    iswitch::ComplexF64 = 1.0 # 1.0 for real time, -im for imaginary time

    # === dimensions and physics
    L::NTuple{D,Float64} # length scales
    N::NTuple{D,Int64}  # grid points in each dimensions
    dV = volume_element(L, N)
    Vol = prod(L)
    ti::Float64 = 0.0    # initial time
    tf::Float64 = 1.0    # final time
    tspan = [ti, tf]
    dt::Float64 = 1e-3 # used for ground state computation

    # === nonlinearity
    g::Float64 = 0.1
    gamma::Float64 = 0.0; @assert gamma >= 0.0 # damping parameter
    mu::Float64 = 0.0 # fixed chemical potential for ground state solution
    sigma2 = init_sigma2(g) 
    # === potential
    params::UserParams = Params() # optional user parameterss
    V0::A = zeros(N)

    # === initial value
    psi_0::A = ones(N) |> complex # initial condition


    # === saving
    nfiles::Bool = false
    Nt::Int64 = 5    # number of saves over (ti,tf)
    t::LinRange{Float64} = LinRange(ti,tf,Nt) # time of saves
    path::String = nfiles ? joinpath(@__DIR__,"data") : @__DIR__
    filename::String = "save"

    # === arrays, transforms, spectral operators
    X::NTuple{D,A} = xvecs(L,N)
    K::NTuple{D,A} = kvecs(L,N)
    T::TransformLibrary{A} = makeT(X,K,A,flags=flags)
    ksquared::A = k2(K, A)
end

InitSim(L,N,A,par) = Sim{length(L), A}(L=L,N=N,params=par)
InitSim(L,N,A) = Sim{length(L), A}(L=L,N=N,params=Params())