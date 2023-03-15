using Pkg
Pkg.activate(".")

using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
using ProgressBars
import JLD2
gr()
GR.usecolorscheme(1)
let sim

L = (50.0,)
N = (256,)

sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

# ====== tiling settings 
tiles = 25

# ====== initialization and unpacking
@unpack_Sim sim
g = -0.587  #corresponds to gamma

# ====== tiling parameters
# in previous simulations:
#   -  max velocity = 1 V_S
#   -  max barrier  = 1.2246 E_S

max_vel = abs(g)     # * 10
max_bar = abs(g)^2 * 1.2246  #* 10

vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
tran = Array{Float64, 2}(undef, (tiles, tiles))
refl = Array{Float64, 2}(undef, (tiles, tiles))


gamma = 0.0
g_param = abs(g) / 2

equation = NPSE
iswitch = 1
x = X[1] |> real
k = K[1] |> real
dV= volume_element(L, N)
reltol = 1e-3
x0 = L[1] / 4
alg = BS3()
maxiters = 50000

@pack_Sim! sim

mask_refl = map(xx -> xx>0, x)
mask_tran = map(xx -> xx<0, x)

iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))

p = plot(x, zeros(length(x)))
for ((vx, vv), (bx, bb)) in ProgressBar(iter)
    @unpack_Sim sim
    @. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)

    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    @. V0 = bb * exp(-(x/0.699)^2)
    if vv == 0.0
        tf = 10.0
    else
        tf = 2*x0/vv
    end
    Nt = 2
    t = LinRange(ti, tf, Nt)
    @info "Computing tile" (vv, bb)
    @pack_Sim! sim

    sol = runsim(sim; info=false)
    #JLD2.@save("tran.jld2", tran)

    if isnothing(sol)
        tran[bx, vx] = NaN
        refl[bx, vx] = NaN
        @info "Computed transmission coefficient" tran[bx, vx]
    else
    final = sol[end]
    # plot!(p, x, abs2.(final))
    
    # time_axis = sol.t
    # u = reduce(hcat, sol.u)
    # u = mapslices(x->xspace(x, sim),u,dims=(1)) 

    # ht = heatmap(real.(x), time_axis, abs2.(u)')
    # display(ht)

    xspace!(final, sim)
    tran[bx, vx] = ns(final, sim, mask_tran)
    refl[bx, vx] = ns(final, sim, mask_refl)
    @info "Computed transmission coefficient" tran[bx, vx]
    end
end

# display(p)

JLD2.@save("tran.jld2", tran)
JLD2.@save("refl.jld2", refl)
norm_bar = bar_list / max_bar
norm_vel = vel_list / max_vel
ht = heatmap(norm_bar, norm_vel, tran')
display(ht)
end

# JLD2.@load "refl.jld2" refl
# norm_bar = bar_list / max_bar
# norm_vel = vel_list / max_vel
# ht = heatmap(norm_bar, norm_vel, tran')
# display(ht)