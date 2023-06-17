using Revise
using Plots, LaTeXStrings
using CondensateDynamics
import Makie, GLMakie
using HDF5, JLD2
using FFTW, CUDA, OrdinaryDiffEq

using LoopVectorization

using ProgressBars, Colors

includet("../examples/plot_axial_evolution.jl")
includet("../examples/plot_isosurfaces.jl")
includet("visual_utils.jl")
includet("simulations_parameters.jl")


```
    file signature:
    tran.JLD2
    refl.JLD2
```
includet("lines.jl")


```
    file signature:
    TODO
```
includet("soliton_ground_state_comparison.jl")


```
    file signature:
    TODO
```
includet("tiles.jl")

static_standard = load_parameters_gs()
dynamic_standard = load_parameters_dy()