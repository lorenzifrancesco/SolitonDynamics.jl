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
# ================ plotting functions
function dense(phi)
    psi = xspace(phi,sim)
    density = abs2.(psi)
    pmax = maximum(density)
    if pmax == 0
        throw("Maximum density is null")
    end
    return density/pmax
end

function isosurface_animation(sol,Nt, sim;file="3Dcollision.gif", framerate=3)
    saveto=joinpath("media",file)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    iter = [Array(xspace(sol[k], sim)) for k in 1:Nt]
    iter = [abs2.(iter[k]) for k in 1:Nt]
    scene = Makie.volume(Makie.@lift(iter[$tindex]/maximum(iter[$tindex])),
                        algorithm =:iso,
                        isovalue=0.2,
                        isorange=0.1,
                        colormap=:greens, 
                        alpha=0.3,
                        transparency=true
    )

    R = 180
    eyeat = Makie.Vec3f0(R,0,0)
    lookat = Makie.Vec3f0(-50,-50,0)
    Makie.record(scene, saveto, 1:Nt; framerate=framerate) do i
        # Makie.update_cam!(scene, eyeat, lookat)
        # Makie.rotate_cam!(scene, 0., -0.4, 0.)
        tindex[] = i
    end
    return
end

function isosurface(sol)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    psol = Array(abs2.(xspace(sol, sim)))
    scene = Makie.volume(psol/maximum(psol),
                        algorithm =:iso,
                        colormap=:greens, 
                        alpha=0.3,
                        isovalue=0.2,
                        isorange=0.1,
    )
    display(scene)
    return
end

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
g = -1.17 * 2*pi
g_param = abs(g)/(4*pi)
reltol = 1e-3
iswitch = 1
vv = 5.0
tf = 3
Nt = 30

x = Array(X[1]) |> real
y = Array(X[2]) |> real
z = Array(X[3]) |> real
dV= volume_element(L, N)

t = LinRange(ti,tf,Nt)
# nfiles = true
maxiters = 5000
x0 = L[1]/4

# tmp = [(exp(-(y^2+z^2)/2) * sqrt(g_param/2)*2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param))) * exp(-im*x*vv) for x in x, y in y, z in z]
tmp = [exp(-(y^2+z^2+(x-x0)^2)/2) * exp(-im*x*vv) for x in x, y in y, z in z]

psi_0 = CuArray(tmp)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()
#1/2*(x^2+y^2+ 3*z^2)
tmp = [ (1/2*(y^2+ z^2) + 0.0*1/2*x^2 + 400*exp(-10*x^2)) for x in x, y in y, z in z]
V0 = CuArray(tmp)
#V(x,y,z,t) = 1/2 * (x^2+y^2+z^2)
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

iter = collect(((collect(enumerate(vel_list[i])), collect(enumerate(bar_list[j]))) for i in 1:tiles for j in 1:tiles))
Threads.@threads for ((vx, vv), (bx, bb)) in ProgressBar(iter)
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
    final = sol[end]
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

end