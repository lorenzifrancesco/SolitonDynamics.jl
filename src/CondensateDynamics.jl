__precompile__()

module CondensateDynamics

# dev settings
using ExportAll

using FFTW, CUDA, Adapt # useful for loading sparse matrix in GPU
using Parameters
using Reexport
using OrdinaryDiffEq, DiffEqCallbacks, SteadyStateDiffEq, DiffEqGPU
using LinearAlgebra, RecursiveArrayTools, LazyArrays
import JLD2 

@reexport using Parameters
import FileIO

export runsim, testsim
export Sim, SISim 
export normalize, printsim

include("methods.jl")
include("types.jl")
include("arrays.jl")
include("solver.jl")
include("normalization.jl")
@exportAll()
end # module CondensateDynamics