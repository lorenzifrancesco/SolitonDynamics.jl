module SolitonDynamics

# dev settings
using ExportAll

using FFTW, CUDA
using Reexport
using OrdinaryDiffEq
import NonlinearSolve, LinearSolve, Interpolations
using LinearAlgebra
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
include("utils.jl")
include("solver.jl")
include("normalization.jl")

@exportAll()

end # module CondensateDynamics
