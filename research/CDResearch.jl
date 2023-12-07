using Revise
using CondensateDynamics
using HDF5, JLD2
using FFTW, CUDA, OrdinaryDiffEq
using Interpolations
using OrderedCollections
using LoopVectorization, LinearAlgebra

# ENV["PYTHON"]="usr/bin/python"
# using PyPlot
using LaTeXStrings
using Plots
import Makie, GLMakie
using ProgressBars, Colors

pyplot()
using PyCall
const plt = pyimport("matplotlib.pyplot")

# # Set other parameters as needed
# plt.rcParams["figure.figsize"] = [8, 6]
# plt.rcParams["font.size"] = 12
# plt.rcParams["lines.linewidth"] = 2
# Set other parameters as needed

includet("_plot_settings.jl")
# pyplot(size=(350, 220))
# if ENV["USER"] == "ubuntu"
#   plotly()
# else
#   pyplot()
# end

includet("../examples/plot_axial_evolution.jl")
includet("../examples/plot_isosurfaces.jl")
includet("visual_utils.jl")
includet("init.jl")
includet("sim_utils.jl")


```
    file signature:
```
includet("solitons.jl")


```
    file signature:
    tran.JLD2
    refl.JLD2
```
includet("lines.jl")


```
    file signature:
```
includet("tiles.jl")

includet("chempot.jl")


includet("aux_collapse.jl")
includet("aux_gs.jl")
includet("aux_collision.jl")
includet("aux_sigma2.jl")
includet("checking.jl")