module SolitonDynamics

# dev settings
using ExportAll

using FFTW, CUDA, Adapt
using Reexport
using OrdinaryDiffEq, DiffEqCallbacks, SteadyStateDiffEq, DiffEqGPU, BoundaryValueDiffEq
import NonlinearSolve, LinearSolve, Interpolations
using LoopVectorization
using LinearAlgebra, RecursiveArrayTools, LazyArrays
using IntervalRootFinding, IntervalArithmetic
import JLD2
using StaticArrays
using ProgressMeter
using Printf

@reexport using Parameters
import FileIO

export runsim, testsim
export Sim, SISim
export normalize, display

include("types.jl")
include("methods.jl")
include("arrays.jl")
include("utils.jl")
include("arrays.jl")
include("solver.jl")
include("normalization.jl")

@exportAll()

end # module CondensateDynamics
