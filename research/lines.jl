using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
using ProgressBars
import JLD2
using Colors

include("simulations_parameters.jl")

plotly(size=(800, 400))

function compute_lines()
    @info "Loading parameters..."

    simulation_dict = get_parameters() 
    include("../src/plot_axial_evolution.jl")
    save_path = "results/"

    file = "tran.pdf"

    @info "Setting phase space..."
    barrier_width = 0.699 # as in SolitonBEC.jl
    max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
    max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY

    lines = 4
    vel_list = LinRange(0.0, max_vel, 4)
    spots = 2
    bar_sweep = LinRange(0.0, max_bar, spots)

    tran = Dict{String, Array{Float64, 2}}
    refl = Dict{String, Array{Float64, 2}}

    p = plot(title = "transmission")

    for (eq_name, eq) in simulation_dict
        @info "Computing lines for $eq_name"
        @info eq_name
        sim = eq
        @unpack_Sim sim

        x = X[1] |> real
        mask_refl = map(xx -> xx>0, x)
        mask_tran = map(xx -> xx<0, x)
        iswitch = 1
        alg = BS3()
        g_param = abs(g) / 2
        x0 = L[1] / 4
        @warn manual
        manual = false

        @pack_Sim! sim
        
        p = plot(x, zeros(length(x)))
        t_eq = Array{Float64, 2}(undef, (lines, spots))
        r_eq = Array{Float64, 2}(undef, (lines, spots))

        for (idv, vel) in enumerate(vel_list) 
            iter = enumerate([[i/spots * bar_sweep[1] + (spots -i)/spots * bar_sweep[2] , 
                            vel] for i in 1:spots])
            avg_iteration_time = 0.0
            
            full_time = @elapsed for (idx, coordinate) in ProgressBar(iter)
                bb = coordinate[1]
                vv = coordinate[2] 
                @info "Computing for $bb and $vv"
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
                @pack_Sim! sim

                avg_iteration_time += @elapsed sol = runsim(sim; info=false)

                if isnothing(sol)
                    t_eq[eq_idv, idx] = NaN
                    r_eq[eq_idv, idx] = NaN
                else
                final = sol.u[end]

                xspace!(final, sim)
                t_eq[idv, idx] = ns(final, sim, mask_tran)
                r_eq[idv, idx] = ns(final, sim, mask_refl)
                end
            end
            @info "Tiling time            = " full_time
            @info "Total time in solver   = " avg_iteration_time
            @info "Average iteration time = " avg_iteration_time / spots^2
            tran[string(eq_name)] = t_eq
            refl[string(eq_name)] = r_eq
            plot!(p, 1:spots, tran[idv, :], label="v = $vel", lw=2, color=:blue)
            display(p)
        end
    end
    return tran, refl
end

tran, refl = compute_lines()
