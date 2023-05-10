__precompile__()

module CondensateDynamics

# dev settings
using ExportAll

using FFTW, CUDA, Adapt # useful for loading sparse matrix in GPU
using Parameters
using Reexport
using OrdinaryDiffEq, DiffEqCallbacks, SteadyStateDiffEq, DiffEqGPU, BoundaryValueDiffEq, NonlinearSolve
using LinearAlgebra, RecursiveArrayTools, LazyArrays
using IntervalRootFinding, IntervalArithmetic
import JLD2 

@reexport using Parameters
import FileIO

export runsim, testsim
export Sim, SISim 
export normalize, printsim

global sigma2_old = Vector{Vector{Float64}}()
global sigma2_new = Vector{Vector{Float64}}()
global time_of_sigma = Vector{Float64}()

include("methods.jl")
include("types.jl")
include("utils.jl")
include("arrays.jl")
include("solver.jl")
include("normalization.jl")
@exportAll()
end # module CondensateDynamics