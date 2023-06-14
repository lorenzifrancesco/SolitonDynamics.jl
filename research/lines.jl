using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
using ProgressBars
import JLD2

plotly(size=(800, 400))

@info "Loading parameters..."
include("simulations_parameters.jl")

simulation_dict = get_parameters() 
sim_gpe_1d    = simulation_dict["GPE_1D_GS"]
sim_npse      = simulation_dict["NPSE_GS"]
sim_gpe_3d    = simulation_dict["GPE_3D_GS"]
sim_npse_plus = simulation_dict["NPSE_plus_GS"]

include("../src/plot_axial_evolution.jl")
save_path = "results/"

file = "tran.pdf"

@info "Setting phase space..."
mask_refl = map(xx -> xx>0, x)
mask_tran = map(xx -> xx<0, x)
p = plot(x, zeros(length(x)))

vel_list = LinRange(0.0, max_vel, 4)
spots = 1e3

for vel in vel_list 
    iter = enumerate([[i/tiles * bar_coord[1] + (tiles -i)/tiles * bar_coord[2] , 
                    i/tiles * vel_coord[1] + (tiles -i)/tiles * vel_coord[2]] for i in 1:tiles])

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


    # JLD2.@load "refl.jld2" refl
    # norm_bar = bar_list / max_bar
    # norm_vel = vel_list / max_vel
    # ht = heatmap(norm_bar, norm_vel, tran')
    # display(ht)
end