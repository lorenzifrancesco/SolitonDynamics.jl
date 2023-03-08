__precompile__()

module CondensateDynamics

# dev settings
using ExportAll

using FFTW, CUDA
using Parameters
using Reexport
using LinearAlgebra
using OrdinaryDiffEq, DiffEqCallbacks, SteadyStateDiffEq
using LazyArrays
import DiffEqGPU
using LinearAlgebra, RecursiveArrayTools

@reexport using Parameters
@reexport using JLD2
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