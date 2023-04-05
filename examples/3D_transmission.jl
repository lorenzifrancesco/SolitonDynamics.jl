using Pkg
Pkg.activate(".")

using CondensateDynamics, OrdinaryDiffEq, LSODA
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
import GR
using CUDA.CUFFT
import Makie, GLMakie

include("plot_isosurfaces.jl")
include("plot_axial_evolution.jl")

# =================== simulation settings
L = (40.0,40.0,40.0)
N = (128,128,128)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

# =================== physical parameters 
@unpack_Sim sim
g = -1.17 * 2 * pi
g_param = abs(g)/(4*pi)
equation = GPE_3D
manual = true
iswitch = 1
x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)
reltol = 1e-3
time_steps = 100
Nt = 30
t = LinRange(ti,tf,Nt)
# nfiles = true
maxiters = 2000
tiles = 8
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY
vx = 4
bx = 4
x0 = L[3]/5
tf = 1
# tf = x0*2/vv

vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
vv = vel_list[vx]
bb = bar_list[bx]

tmp = [exp(-(x^2+y^2)/2) * sqrt(g_param/2) * 2/(exp(g_param*(z-x0)) + exp(-(z-x0)*g_param)) * exp(-im*(z-x0)*vv) for x in x, y in y, z in z]
psi_0 = CuArray(tmp)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()
tmp = [1/2*(x^2+y^2) + bb*exp(-(x/barrier_width)^2) for x in x, y in y, z in z]
V0 = CuArray(tmp)
#V(x,y,z,t) = 1/2 * (x^2+y^2+z^2)
@pack_Sim! sim


# ===================== simulation
sol = runsim(sim; info=false)
u = sol.u
# =================== plotting and collect 

plot_axial_heatmap(u, sol.t, sim, 3)
plot_final_density(u, psi_0, sim, 3; info=true)

@info "Building animation..."
isosurface_animation(sol,length(sol), sim; framerate=5)
@info "Completed."