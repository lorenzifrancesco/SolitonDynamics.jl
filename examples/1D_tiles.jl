@time begin # 30s startup
using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
using ProgressBars
import JLD2
end

gr()
GR.usecolorscheme(1)

file = "tran.pdf"
let sim
saveto=joinpath("media/1D",file)

L = (50.0,)
N = (256,)

sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

# ====== tiling parameters
# in previous simulations:
#   -  max velocity = 1 V_S
#   -  max barrier  = 1.2246 E_S
# In SolitonBEC: g = -1.17

# ====== initialization and unpacking
@unpack_Sim sim
# ======= simulation custom parameters
equation = NPSE
solver = SplitStep 
g = -1.17  #corresponds to gamma -0.587
gamma = 0.0
tiles = 3
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY
reltol = 1e-4
abstol = 1e-4
x0 = L[1] / 4

# other computations
vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
tran = Array{Float64, 2}(undef, (tiles, tiles))
refl = Array{Float64, 2}(undef, (tiles, tiles))

iswitch = 1
g_param = abs(g) / 2

x = X[1] |> real
k = K[1] |> real
dV= volume_element(L, N)
sigma2 = init_sigma2(g)

alg = BS3()
maxiters = 50000

@pack_Sim! sim

mask_refl = map(xx -> xx>0, x)
mask_tran = map(xx -> xx<0, x)
p = plot(x, zeros(length(x)))

# iter = collect(((collect(enumerate(vel_list))[i], collect(enumerate(bar_list))[j]) for i in 1:tiles for j in 1:tiles))
# display(iter)
# @assert Threads.nthreads() > 1
# Threads.@threads for ((vx, vv), (bx, bb)) in iter
iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
for ((vx, vv), (bx, bb)) in ProgressBar(iter)
    @unpack_Sim sim
    @. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)

    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    @. V0 = bb * exp(-(x/barrier_width)^2 /2)
    if vv == 0.0
        tf = 10.0
    else
        tf = 2*x0/vv
    end
    Nt = 2
    t = LinRange(ti, tf, Nt)
    #@info "Computing tile" (vv, bb)
    @pack_Sim! sim

    @time sol = runsim(sim; info=false)

    if isnothing(sol)
        tran[bx, vx] = NaN
        refl[bx, vx] = NaN
        @info "T = " tran[bx, vx]
    else
    #JLD2.@save("tran.jld2", tran)
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
    @info "T = " tran[bx, vx]
    #@info "difference wrt NPSE alone: " tran[bx, vx] - mat[bx, vx]
    print("\n")
    end

end

# display(p)

# JLD2.@save("tran.jld2", tran)
# JLD2.@save("refl.jld2", refl)
norm_bar = bar_list / max_bar
norm_vel = vel_list / max_vel
ht = heatmap(norm_bar, norm_vel, tran')
display(ht)
savefig(ht, saveto)
end

# JLD2.@load "refl.jld2" refl
# norm_bar = bar_list / max_bar
# norm_vel = vel_list / max_vel
# ht = heatmap(norm_bar, norm_vel, tran')
# display(ht)