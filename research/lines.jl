function all_lines(precomputed=false)
    plotly(size=(800, 400))
    saveto = "media/lines.pdf"
    # Computing all the lines
    # if files are already saved, load them

    if precomputed && isfile("lines_tran.jld2") && isfile("lines_refl.jld2")
        tran = JLD2.load("lines_tran.jld2", "tran")
        refl = JLD2.load("lines_refl.jld2", "refl")
    else
        tran, refl = compute_lines()
        JLD2.save("lines_tran.jld2", "tran", tran)
        JLD2.save("lines_refl.jld2", "refl", refl)
    end

    # Plotting
    p = plot(title = "transmission")
    for (eq, data) in tran
        for i in eachindex(data[:, 1])
            plot!(p, 1:length(data[i, :]), data[i, :], label=eq, lw=2)
        end
    end
    display(p)
    savefig(p, saveto)
    return nothing
end

function compute_lines(max_vel=0.1, max_bar=1.68, lines=4, spots=10)
    @info "Loading parameters..."
    simulation_dict = load_parameters_dy() # FIXME no import inside functions
    @info "Setting phase space..."
    barrier_width = 0.699 # as in SolitonBEC.jl
    
    vel_list = LinRange(0.0, max_vel, 4)
    bar_sweep = LinRange(0.0, max_bar, spots)
    
    tran = Dict{String, Array{Float64, 2}}()
    refl = Dict{String, Array{Float64, 2}}()
    
    p = plot(title = "transmission")
    
    for (eq_name, eq) in ProgressBar(simulation_dict)
        @info "Computing lines for $eq_name"
        @info eq_name
        sim_0 = eq
        @unpack_Sim sim_0
        
        x = X[1] |> real
        mask_refl = map(xx -> xx>0, x)
        mask_tran = map(xx -> xx<0, x)
        iswitch = 1
        alg = BS3()
        g_param = abs(g) / 2
        x0 = L[1] / 4
        manual = false
        
        @pack_Sim! sim_0
        
        p = plot(x, zeros(length(x)))
        t_eq = Array{Float64, 2}(undef, (lines, spots))
        r_eq = Array{Float64, 2}(undef, (lines, spots))
        
        for (idv, vel) in ProgressBar(enumerate(vel_list))
            iter = enumerate([[i/spots * bar_sweep[1] + (spots -i)/spots * bar_sweep[2] , 
            vel] for i in 1:spots])
            avg_iteration_time = 0.0
            full_time = @elapsed for (idx, coordinate) in collect(ProgressBar(iter))
                bb = coordinate[1]
                vv = coordinate[2] 
                # @info "Computing for $bb and $vv"
                sim = deepcopy(sim_0)
                
                # initialize sim
                @unpack_Sim sim
                
                psi_0 .= psi_0 * 0
                @. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)
                psi_0 .= psi_0 / sqrt(ns(psi_0, sim))
                
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
            plot!(p, 1:spots, t_eq[idv, :], label="v = $vel", lw=2, color=:blue)
            display(p)
        end
    end
    return tran, refl
end