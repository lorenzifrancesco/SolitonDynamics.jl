module SolitonDynamics

# library and dev utils  
using ExportAll, Reexport

# Numerics
using FFTW, CUDA
import NonlinearSolve, LinearSolve, Interpolations
using LinearAlgebra
using IntervalRootFinding, IntervalArithmetic
using StaticArrays

# IO
import JLD2, FileIO, CSV
using ProgressMeter
using Printf
@reexport using Parameters

export runsim, testsim
export Sim, SISim
export normalize, display

include("methods.jl")
include("types.jl")
include("utils.jl")
include("analytics.jl")
include("functionals.jl")
include("sigma.jl")
include("solver.jl")
include("normalization.jl")
include("save.jl")

@exportAll()

end # module CondensateDynamics
