using Pkg
Pkg.activate(".")

using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
using ProgressBars
gr()
GR.usecolorscheme(1)

L = (40.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

# ====== tiling settings 
tiles = 50
vel_list = LinRange(0, 0.3, tiles)
bar_list = LinRange(0, 1, tiles)
tran = Array{Float64, 2}(undef, (tiles, tiles))
refl = Array{Float64, 2}(undef, (tiles, tiles))

# ====== initialization and unpacking
@unpack_Sim sim
g = -3.0
equation = GPE_1D
iswitch = 1
x = X[1] |> real
k = K[1] |> real
dV= volume_element(L, N)
reltol = 1e-4
tf = 10.0
x0 = -10
alg = Tsit5()
@pack_Sim! sim

mask_refl = map(xx -> xx>0, x)
mask_tran = map(xx -> xx<0, x)

iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))

for ((vx, vv), (bx, bb)) in ProgressBar(iter)
    @unpack_Sim sim
    @. psi_0 = exp(-(x-x0)^2/2) * exp(im*x*vv)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    @. V0 = bb * exp(-10*x^2)
    @pack_Sim! sim

    sol = runsim(sim; info=false)

    final = sol[end]
    xspace!(final, sim)
    tran[bx, vx] = ns(final, sim, mask_tran)
    refl[bx, vx] = ns(final, sim, mask_refl)
end
@save("tran.jld2", tran)
@save("refl.jld2", refl)

ht = heatmap(1:tiles, 1:tiles, tran)
display(ht)