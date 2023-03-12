using Pkg
Pkg.activate(".")

using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
gr()
GR.usecolorscheme(1)

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (50.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
g = -0.587
gamma = 0.0
Nt = 10
t = LinRange(ti,tf,Nt)
g_param = abs(g) / 2
equation = GPE_1D
iswitch = 1
x = X[1] |> real
k = K[1] |> real
dV= volume_element(L, N)
reltol = 1e-9
x0 = L[1] / 3

# ====== tiling settings 
tiles = 20
vel_list = LinRange(0, 1.17, tiles)
bar_list = LinRange(0, 1.685, tiles)

vv = vel_list[1]
bb = bar_list[1]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
maxiters = 50000
#@. psi_0 = exp(-x^2/2) * exp(-im*x*10)
@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*x*vv)

psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS5()
@. V0 = bb*exp(-(x/0.699)^2)

@pack_Sim! sim

mask_tran = map(xx -> xx<0, x)


sol = runsim(sim; info=false)
time_axis = sol.t
u = reduce(hcat, sol.u)
@info "time steps" size(time_axis)
final = u[:, end]
xspace!(final, sim)
xspace!(psi_0, sim)
@info "final distribution norm squared: " ns(final, sim)

p = plot(real.(x), abs2.(psi_0), label="initial")
plot!(p, real.(x), abs2.(final), label="final")

@info "transmitted norm" ns(final, sim, mask_tran)

u = mapslices(x->xspace(x, sim),u,dims=(1)) 

#map(x -> xspace!(x, sim), sol)
ht = heatmap(real.(x), time_axis, abs2.(u)')
display(ht)
