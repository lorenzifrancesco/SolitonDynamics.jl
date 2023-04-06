using Pkg
Pkg.activate(".")

using CondensateDynamics, OrdinaryDiffEq, LSODA
import JLD2
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
using CUDA.CUFFT
import Makie, GLMakie
using ProgressBars

JULIA_CUDA_SOFT_MEMORY_LIMIT ="95%"
include("plot_isosurfaces.jl")
include("plot_axial_evolution.jl")

file = "3Dtran.pdf"
let sim
saveto=joinpath("media/1D",file)

# =================== simulation settings
L = (40.0,40.0,40.0)
N = (128,128,128)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

# =================== physical parameters 
@unpack_Sim sim
equation = GPE_3D
manual = true
g = -1.17 # this is the 1D g
g_param = abs(g)/2
reltol = 1e-3
iswitch = 1
vv = 5.0
tf = 3
Nt = 30
dt = 0.1

x = Array(X[1]) |> real
y = Array(X[2]) |> real
z = Array(X[3]) |> real
dV= volume_element(L, N)

t = LinRange(ti,tf,Nt)
# nfiles = true
maxiters = 10000
x0 = L[1]/4

kspace!(psi_0, sim)
alg = BS3()

@pack_Sim! sim


# ===================== tiling
tiles = 4
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY

vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
tran = Array{Float64, 2}(undef, (tiles, tiles))
refl = Array{Float64, 2}(undef, (tiles, tiles))

mask_refl = map(xx -> xx>0, CuArray(x))
mask_tran = map(xx -> xx<0, CuArray(x))

# iter = collect(((collect(enumerate(vel_list[i])), collect(enumerate(bar_list[j]))) for i in 1:tiles for j in 1:tiles))
# Threads.@threads for ((vx, vv), (bx, bb)) in ProgressBar(iter)
iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
for ((vx, vv), (bx, bb)) in ProgressBar(iter)
    @info "Computing tile" (vv, bb)

    # ===================== tile simulation parameters
    @unpack_Sim sim
    tmp = [ (1/2*(y^2+ z^2) + 0.0*1/2*x^2 + bb*exp(-(x/barrier_width)^2)) for x in x, y in y, z in z]
    V0 = CuArray(tmp)

    if vv == 0.0
        tf = 10.0
    else
        tf = 2*x0/vv
    end
    Nt = 2
    t = LinRange(ti, tf, Nt)

    tmp = [exp(-(y^2+z^2+(x-x0)^2)/2) * exp(-im*x*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    initial_state = psi_0
    kspace!(psi_0, sim)
    @pack_Sim! sim
    
    @info "Running solver..."
    sol = runsim(sim; info=false)
    final = sol.u[end]
    @info "Run complete, computing transmission..."
    xspace!(final, sim)
    tran[bx, vx] = ns(final, sim, mask_tran)
    refl[bx, vx] = ns(final, sim, mask_refl)
    @info "T = " tran[bx, vx]
end

JLD2.@save("tran.jld2", tran)
JLD2.@save("refl.jld2", refl)
norm_bar = bar_list / max_bar
norm_vel = vel_list / max_vel
ht = heatmap(norm_bar, norm_vel, tran')
display(ht)
savefig(ht, saveto)

display(sol.destats)
display(sol.retcode)
end