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


L = (70.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)


# ======= packing sim
@unpack_Sim sim
g = -0.0
gamma = 0.0

g_param = abs(g) /2
equation = GPE_1D
solver = SplitStep 
sigma2 = init_sigma2(g)

iswitch = 1
x = X[1] |> real
k = K[1] |> real
dV= volume_element(L, N)
reltol = 1e-4
abstol = 1e-4
x0 = L[1] / 5

# ====== tiling settings 
tiles = 25
vel_list = LinRange(0, abs(g), tiles)
bar_list = LinRange(0, abs(g)^2, tiles)

vv = vel_list[15]
bb = bar_list[15]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
Nt = 200
t = LinRange(ti, tf, Nt)


time_steps = 200
dt = (tf-ti)/time_steps
maxiters = time_steps * 2

#@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)
@. psi_0 = exp(-(x-x0)^2/2)*exp(-im*30*x)
psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()
#@. V0 = bb*exp(-(x/0.699)^2)
@. V0 = 1/2*x^2
@pack_Sim! sim


mask_tran = map(xx -> xx<0, x)

sol = runsim(sim; info=false)
if isnothing(sol)
    throw("NPSE collapse detected, cannot proceed further to plots...")
end
time_axis = sol.t |> real
u = reduce(hcat, sol.u)

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
@info "max" g*maximum(abs2.(u))

display(sol.destats)
display(sol.retcode)