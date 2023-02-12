module CondensateDynamics
import FFTW, CUDA
using Parameters

include("solver.jl")
include("types.jl")
include("arrays.jl")
end # module CondensateDynamics