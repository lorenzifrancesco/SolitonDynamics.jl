using CondensateDynamics, OrdinaryDiffEq, LSODA
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
import GR
using CUDA.CUFFT
import Makie, GLMakie
import JLD2

include("plot_axial_evolution.jl")
include("plot_isosurfaces.jl")
save_path = "results/"
use_precomputed = false
# ================ plotting functions

gr()
GR.usecolorscheme(1)

# =================== simulation settings
L = (40.0,20.0,20.0)
N = (512, 64, 64)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

# =================== physical parameters
@unpack_Sim sim

# "collapse is visible "

g_param = 0.4
g = - g_param * 4 * pi

gamma = 0.0
mu = 0.0

equation = GPE_3D
solver = SplitStep
manual = true
iswitch = -im
reltol = 1e-5
x0 = 0.0
vv = 0.0


x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)
tf = 2.0

Nt = 30
t = LinRange(ti,tf,Nt)
maxiters = 0
abstol = 1e-6
dt = 0.02


tmp = [exp(-((x-x0)^2+y^2+z^2)/2) * exp(-im*x*vv) for x in x, y in y, z in z]
psi_0 = CuArray(tmp)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()

tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
V0 = CuArray(tmp)
@pack_Sim! sim

# ===================== simulation
@info "computing GPE_3D" 
if isfile(join([save_path, "3d_gs.jld2"])) && use_precomputed
    @info "\t using precomputed solution 3d_gs.jld2" 
    JLD2.@load join([save_path, "3d_gs.jld2"]) u
else
    sol = runsim(sim; info=true)
    u = sol.u
    JLD2.@save join([save_path, "3d_gs.jld2"]) u
end

# =================== plotting and collect 
plot_final_density([u], sim, 1; info=true, label="final")


# transverse view
aa = Array(abs2.(xspace(initial_state, sim)))
axial_density = sum(aa, dims=(2, 3))[:, 1, 1] * sim.dV
q = heatmap(aa[190, :, :] / axial_density[190])
display(q)

# top view
# aa = Array(abs2.(xspace(initial_state, sim)))
# qt = heatmap(aa[:, :, 32])
# display(qt)

# s2 = estimate_sigma2(u, sim)
# @info "estimated sigma2" s2
# @info "minimum   sigma2" minimum(s2)
# q = plot(sim.X[1] |> real, s2, label="bella zio")
# display(q)
#@info "Building animation..."
#isosurface_animation(u, length(u), sim; framerate=5)
@info "Completed."