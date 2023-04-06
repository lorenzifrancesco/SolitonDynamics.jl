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

# ====== initialization and unpacking
@unpack_Sim sim
# ======= simulation custom parameters
equation = GPE_1D
solver = SplitStep 
manual = false
time_steps = 2000
g = -1.17  #corresponds to gamma -0.587
gamma = 0.0
tiles = 100
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY
reltol = 1e-4
abstol = 1e-4
x0 = L[1] / 4

# other computations
vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
tran = Array{Float64, 1}(undef, (tiles))
refl = Array{Float64, 1}(undef, (tiles))


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

vel_coord = [0 , 1] * max_vel

bar_coord = [0.2 , 0.2] * max_bar 

iter = enumerate([[i/tiles * bar_coord[1]- (tiles -i)/tiles * bar_coord[2] , 
                   i/tiles * vel_coord[1]- (tiles -i)/tiles * vel_coord[2]] for i in 1:tiles])
avg_iteration_time = 0.0
full_time = @elapsed for (idx, coordinate) in ProgressBar(iter)
    bb = coordinate[1]
    vv = coordinate[2] 
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
    dt = (tf-ti)/time_steps
    #@info "Computing tile" (vv, bb)
    @pack_Sim! sim

    avg_iteration_time += @elapsed sol = runsim(sim; info=false)

    if isnothing(sol)
        tran[idx] = NaN
        refl[idx] = NaN
        #@info "T = " tran[bx, vx]
    else
    #JLD2.@save("tran.jld2", tran)
    final = sol.u[end]

    xspace!(final, sim)
    tran[idx] = ns(final, sim, mask_tran)
    refl[idx] = ns(final, sim, mask_refl)
    #@info "T = " tran[bx, vx]
    #@info "difference wrt NPSE alone: " tran[bx, vx] - mat[bx, vx]
    print("\n")
    end
end
@info "Tiling time            = " full_time
@info "Total time in solver   = " avg_iteration_time
@info "Average iteration time = " avg_iteration_time / tiles^2

# display(p)

# JLD2.@save("tran.jld2", tran)
# JLD2.@save("refl.jld2", refl)
norm_bar = bar_list / max_bar
norm_vel = vel_list / max_vel

p = plot(1:tiles, tran, label="transmission")
display(p)
end

# JLD2.@load "refl.jld2" refl
# norm_bar = bar_list / max_bar
# norm_vel = vel_list / max_vel
# ht = heatmap(norm_bar, norm_vel, tran')
# display(ht)