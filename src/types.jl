abstract type TransformLibrary{A<:AbstractArray} end
abstract type Space end
# abstract parameter type: can also be a function of simulation time

abstract type PotentialType end
abstract type Gaussian <: PotentialType end
abstract type Elliptical <: PotentialType end

struct EquationType
  name::String
  number::Int64
  D::Int64
  color::Symbol
  linestyle::Symbol
end

const GPE_1D = EquationType("G1", 1, 1, :grey, :dashdot)
const NPSE = EquationType("N", 2, 1, :green, :dot)
const NPSE_plus = EquationType("Np", 3, 1, :green, :dash)
const GPE_3D = EquationType("G3", 4, 3, :red, :solid)
const CQGPE = EquationType("CQ", 5, 1, :grey, :dot)

function isless(eq1::EquationType, eq2::EquationType)
  return eq1.name < eq2.name
end

abstract type Simulation{D} end

struct Solver
  number::Int64
  spectral::Bool
end

const SplitStep = Solver(1, true)
const CrankNicholson = Solver(2, false)
const PredictorCorrector = Solver(3, false)
const BackwardEuler = Solver(4, false)


abstract type UserParams end
abstract type Method end

struct NpseCollapse <: Exception
  var::Float64
end
struct Gpe3DCollapse <: Exception
  var::Float64
end

Base.showerror(io::IO, e::NpseCollapse) =
  print(io, "NPSE collapse detected, g * max(|f|^2) = ", e.var, "!")

@with_kw mutable struct Params <: UserParams
  κ = 0.0 # a placeholder
end

struct Transforms{T} <: TransformLibrary{T}
  Txk::AbstractFFTs.ScaledPlan{ComplexF64,FFTW.cFFTWPlan{ComplexF64,-1,false,1,UnitRange{Int64}},Float64}
  Txk!::AbstractFFTs.ScaledPlan{ComplexF64,FFTW.cFFTWPlan{ComplexF64,-1,true,1,UnitRange{Int64}},Float64}
  Tkx::AbstractFFTs.ScaledPlan{ComplexF64,FFTW.cFFTWPlan{ComplexF64,1,false,1,UnitRange{Int64}},Float64}
  Tkx!::AbstractFFTs.ScaledPlan{ComplexF64,FFTW.cFFTWPlan{ComplexF64,1,true,1,UnitRange{Int64}},Float64}
end

@with_kw mutable struct GPUTransforms{A} <: TransformLibrary{A}
  Txk::AbstractFFTs.ScaledPlan{ComplexF64,CUDA.CUFFT.cCuFFTPlan{ComplexF64,-1,false,3},Float64}
  Txk!::AbstractFFTs.ScaledPlan{ComplexF64,CUDA.CUFFT.cCuFFTPlan{ComplexF64,-1,true,3},Float64}
  Tkx::AbstractFFTs.ScaledPlan{ComplexF64,CUDA.CUFFT.cCuFFTPlan{ComplexF64,1,false,3},Float64}
  Tkx!::AbstractFFTs.ScaledPlan{ComplexF64,CUDA.CUFFT.cCuFFTPlan{ComplexF64,1,true,3},Float64}
  #psi::ArrayPartition = crandnpartition(D,N,A)
end

@with_kw mutable struct Sim{D,A<:AbstractArray}
  name::String = "default"
  # === solver and algorithm
  equation::EquationType = GPE_1D
  manual::Bool = false
  solver::Solver = SplitStep
  graphics::Bool = false
  alg = Tsit5() # default solver
  reltol::Float64 = 1e-3  # default tolerance; may need to use 1e-7 for corner cases
  abstol::Float64 = 1e-3
  maxiters::Int64 = 5000
  flags::UInt32 = FFTW.MEASURE # choose a plan. PATIENT, NO_TIMELIMIT, EXHAUSTIVE
  iswitch::ComplexF64 = 1.0 # 1.0 for real time, -im for imaginary time

  # === dimensions and physics
  L::NTuple{D,Float64} # length scales
  N::NTuple{D,Int64}  # grid points in each dimensions
  dV::Float64 = volume_element(L, N)
  Vol::Float64 = prod(L)
  ti::Float64 = 0.0    # initial time
  tf::Float64 = 1.0    # final time
  tspan = [ti, tf]
  dt::Float64 = 1e-3  # used in manual solvers
  time_steps = 5000   # used in manual solvers
  # === nonlinearity
  g::Float64 = 0.0
  gamma_damp::Float64 = 0.0
  @assert gamma_damp >= 0.0 # damping parameter
  mu::Float64 = 0.0 # fixed chemical potential for ground state solution
  sigma2 = init_sigma2(g)
  collapse_threshold::Float64 = 0.1
  # === potential
  params::UserParams = Params() # optional user parameterss
  V0::A = zeros(N)

  # === initial value
  psi_0::A = ones(N) |> complex # initial condition

  # === saving
  nfiles::Bool = false
  Nt::Int64 = 5    # number of saves over (ti,tf)
  t::LinRange{Float64} = LinRange(ti, tf, Nt) # time of saves
  path::String = nfiles ? joinpath(@__DIR__, "data") : @__DIR__
  filename::String = "save"

  # === arrays, transforms, spectral operators
  X::Vector{A} = xvecs(L, N)
  K::Vector{A} = kvecs(L, N)
  T::TransformLibrary{A} = makeT(X, K, A, flags=flags)
  ksquared::A = k2(K, A)

  # === graphics
  color::Symbol = get_color(equation)
  linestyle::Symbol = get_linestyle(equation)
end

function isless(sim1::Sim, sim2::Sim)
  return isless(sim1.equation, sim2.equation)
end

@with_kw mutable struct CustomSolution
  u::AbstractArray
  sigma::AbstractArray=[0.0]
  t::Any
  cnt::Int64 = 0
end

function init_sim(L, N)
  if length(L) == 1
    sim = Sim{1,Array{ComplexF64}}(L=L, N=N)
  elseif length(L) == 3
    sim = Sim{3,CuArray{ComplexF64}}(L=L, N=N)
  else
    throw(MethodError)
  end
  sim
end

import Base.display
function display(s::Sim)
  @printf("Simulation in D = %i, equation = %s", length(s.N), s.equation.name)
  if length(s.N) == 1
    @printf("\n- Domain: \n\tL = %5.1f, \n\tN = %5i", s.L[1], s.N[1])
  else
    @printf("\n- Domain: \n\tL = (%5.1f, %5.1f, %5.1f), \n\tN = (%5i, %5i, %5i)", s.L[1], s.L[2], s.L[3], s.N[1], s.N[2], s.N[3])
  end
  @printf("\n- Nonlinearity: \n\t           g = %4.3f, \n\tcollapse_thr = %5.3f", g2gamma(s.g, s.equation), s.collapse_threshold)
  @printf("\n- %s time", s.iswitch == 1 ? "Real" : "Imag")
  @printf("\n- dt = %4.3f, maxiters = %4i", s.dt, s.maxiters)
  @printf("\n- Number of saves = %3i", s.Nt)
  nothing
end

function get_color(eq::EquationType)
  if eq == GPE_3D
    return :red
  elseif eq == GPE_1D
    return :grey
  elseif eq == NPSE
    return :green
  elseif eq == NPSE_plus
    return :green
  elseif eq == CQGPE
    return :blue
  else
    return :black
  end
end

function get_linestyle(eq::EquationType)
  if eq == GPE_3D
    return :solid
  elseif eq == GPE_1D
    return :dashdot
  elseif eq == NPSE
    return :dot
  elseif eq == NPSE_plus
    return :dash
  elseif eq == CQGPE
    return :dashdot
  else
    return :solid
  end
end

function paper_name(eq::EquationType)
  if eq == GPE_3D
    return "3D-GPE"
  elseif eq == GPE_1D
    return "1D-GPE"
  elseif eq == NPSE
    return "NPSE"
  elseif eq == NPSE_plus
    return "NPSE+"
  elseif eq == CQGPE
    return "CQ-GPE"
  else
    throw("equation not found")
  end
end