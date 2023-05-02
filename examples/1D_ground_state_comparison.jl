using Pkg
Pkg.activate(".")

using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
using LSODA
import CondensateDynamics.V
import FFTW
import JLD2
gr()
GR.usecolorscheme(1)
include("plot_axial_evolution.jl")
save_path = "results/"
# =========================================================
## 1D-GPE 

L = (40.0,)
N = (1024,)
sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim_gpe_1d
iswitch = -im
equation = GPE_1D
manual = true
solver = SplitStep
time_steps = 0

# interaction parameter
g_param = 0.0
use_precomputed = false
# ============
maxiters = 100000

g = - 2 * g_param

n = 100
as = g_param / n
abstol = 1e-6

dt = 0.022

x = X[1]
k = K[1]
dV= volume_element(L, N)

flags = FFTW.EXHAUSTIVE

width = 7
tf = Inf
# SPR condensate bright soliton t in units of omega_perp^-1
analytical_gs = zeros(N)
@. analytical_gs = sqrt(g_param/2) * 2/(exp(g_param*x) + exp(-x*g_param))

psi_0 .= exp.(-(x/10).^2)
psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
initial_state .= psi_0

kspace!(psi_0, sim_gpe_1d)
#@. V0 = 1/2*0.01*x^2

@pack_Sim! sim_gpe_1d

# =========================================================
## NPSE 

sim_npse = sim_gpe_1d
@unpack_Sim sim_npse
equation = NPSE
@info "gamma: " g_param
if g_param > 2/3
    @warn "we should expect NPSE collapse"
end
sigma2 = init_sigma2(g)
@pack_Sim! sim_npse

# =========================================================
## 3D-GPE 

L = (40.0,40.0,40.0)
N = (128, 128, 128)
sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim_gpe_3d
iswitch = -im
equation = GPE_3D
manual = true
solver = SplitStep

time_steps = 10000
g = - g_param * (4*pi)

abstol = 1e-6
maxiters = 1000
dt = 0.022
x0 = 0.0
vv = 0.0

x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)

flags = FFTW.EXHAUSTIVE

tf = Inf
tmp = [exp(-((x-x0)^2+y^2+z^2)/2) * exp(-im*x*vv) for x in x, y in y, z in z]
psi_0 = CuArray(tmp)
psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))

kspace!(psi_0, sim_gpe_3d)
tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
V0 = CuArray(tmp)

@pack_Sim! sim_gpe_3d
# =========================================================
@info "computing GPE_1D" 
if isfile(join([save_path, "gpe_1d.jld2"])) && use_precomputed
    @info "\t using precomputed solution gpe_1d.jld2" 
    JLD2.@load join([save_path, "gpe_1d.jld2"]) gpe_1d
else
    sol = runsim(sim_gpe_1d; info=false)
    gpe_1d = sol.u
    JLD2.@save join([save_path, "gpe_1d.jld2"]) gpe_1d
end

@info "computing NPSE" 
if isfile(join([save_path, "npse.jld2"])) && use_precomputed
    @info "\t using precomputed solution npse.jld2" 
    JLD2.@load join([save_path, "npse.jld2"]) npse
else
    sol = runsim(sim_npse; info=false)
    npse = sol.u
    JLD2.@save join([save_path, "npse.jld2"]) npse
end

# @info "computing GPE_3D" 
# if isfile(join([save_path, "gpe_3d.jld2"])) && use_precomputed
#     @info "\t using precomputed solution gpe_3d.jld2" 
#     JLD2.@load join([save_path, "gpe_3d.jld2"]) gpe_3d
# else
#     sol = runsim(sim_gpe_3d; info=false)
#     gpe_3d = sol.u
#     JLD2.@save join([save_path, "gpe_3d.jld2"]) gpe_3d
# end

p = plot_final_density([gpe_1d], sim_gpe_1d; label="GPE_1D")
#plot_final_density!(p, [analytical_gs], sim_gpe_1d; label="analytical_GPE_1D", doifft=false)
plot_final_density!(p, [npse], sim_npse; label="NPSE")

# plot_final_density!(p, [gpe_3d], sim_gpe_3d, 3; label="GPE_3D")

display(p)