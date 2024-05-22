module SolitonDynamics

# library and dev utils  
using ExportAll, Reexport

# Numerics
using FFTW, CUDA
using OrdinaryDiffEq
import NonlinearSolve, LinearSolve, Interpolations
using LinearAlgebra
using IntervalRootFinding, IntervalArithmetic
using StaticArrays

# IO
import JLD2, FileIO
using ProgressMeter
using Printf
@reexport using Parameters

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
