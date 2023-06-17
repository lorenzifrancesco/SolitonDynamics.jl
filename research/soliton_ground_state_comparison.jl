function get_ground_states()
    @assert false
    # TODO plotly(size=(800, 400))
    @info "Loading parameters..."
    simulation_dict = get_parameters() 

    save_path = "results/"
    gamma_param_list = [0.15, 0.4, 0.6]
    use_precomputed = false
    take_advantage = true
    @info "Starting simulations..."

    for gamma_param in gamma_param_list
        # update simulation parameters
        sim_gpe_1d.g    = -2 * gamma_param 
        sim_npse.g      = -2 * gamma_param
        sim_npse_plus.g = -2 * gamma_param
        sim_gpe_3d.g    = - 4 * pi * gamma_param
        
        # =========================================================
        Plots.CURRENT_PLOT.nullableplot = nothing
        p = plot_final_density([analytical_gs], sim_gpe_1d; label="analytical", color=:orange, doifft=false, ls=:dashdot, title="gamma = $gamma_param")
        plot_final_density!(p, [initial_state_gpe_1d], sim_gpe_1d; label="initial_GPE_1D", color=:grey, doifft=false, ls=:dashdot)
        
        add_string = [join(["_gamma_", gamma_param])]

        @info "computing GPE_1D" 
        if isfile(join([save_path, join(["gpe_1d",add_string, ".jld2"])])) && use_precomputed
            @info "\t using precomputed solution gpe_1d.jld2" 
            gpe_1d = JLD2.load(join([save_path, join(["gpe_1d",add_string, ".jld2"])]))["gpe_1d"]
        else
            sol = runsim(sim_gpe_1d; info=true)
            gpe_1d = sol.u
            JLD2.save(join([save_path, join(["gpe_1d",add_string, ".jld2"])]), "gpe_1d",  gpe_1d)
        end
        plot_final_density!(p, [gpe_1d], sim_gpe_1d; label="GPE_1D", color=:blue, ls=:dash)

        # estimate width
        initial_sigma_improved = sqrt(sum(abs2.(sim_gpe_1d.X[1])  .* abs2.(gpe_1d) * dV) |> real)
        @warn "initial sigma improved" initial_sigma_improved
        if take_advantage 
            sim_npse.psi_0 = gpe_1d
        end
        @info "computing NPSE"
        if isfile(join([save_path, join(["npse",add_string, ".jld2"])])) && use_precomputed
            @info "\t using precomputed solution npse.jld2" 
            npse = JLD2.load(join([save_path, join(["npse",add_string, ".jld2"])]))["npse"]
        else
            sol = runsim(sim_npse; info=true)
            npse = sol.u
            JLD2.save(join([save_path, join(["npse",add_string, ".jld2"])]), "npse",  npse)
        end
        plot_final_density!(p, [npse], sim_npse; label="NPSE", color=:green, ls=:dot)

        if take_advantage
            sim_npse_plus.psi_0 = npse
        end
        @info "computing NPSE_plus" 
        if isfile(join([save_path, join(["npse_plus",add_string, ".jld2"])])) && use_precomputed
            @info "\t using precomputed solution npse_plus.jld2" 
            npse_plus = JLD2.load(join([save_path, join(["npse_plus",add_string, ".jld2"])]))["npse_plus"]
        else
            sol = runsim(sim_npse_plus; info=true)
            npse_plus = sol.u
            JLD2.save(join([save_path, join(["npse_plus",add_string, ".jld2"])]), "npse_plus",  npse_plus)
        end
        plot_final_density!(p, [npse_plus], sim_npse_plus; label="NPSE_der", ls=:dash, color=:green)


        if take_advantage
            x = Array(sim_gpe_3d.X[1])
            y = Array(sim_gpe_3d.X[2])
            z = Array(sim_gpe_3d.X[3])
            initial_sigma_improved *= 5
            tmp = [exp(-(((x-x0)/initial_sigma_improved)^2 + (y^2 + z^2)/2)) * exp(-im*x*vv) for x in x, y in y, z in z]
            sim_gpe_3d.psi_0 = CuArray(tmp)
            sim_gpe_3d.psi_0 .= sim_gpe_3d.psi_0 / sqrt(sum(abs2.(sim_gpe_3d.psi_0) * sim_gpe_3d.dV))
            initial_3d = copy(sim_gpe_3d.psi_0)
            kspace!(sim_gpe_3d.psi_0, sim_gpe_3d)
        end
        @info "computing GPE_3D" 
        if isfile(join([save_path, join(["gpe_3d",add_string, ".jld2"])])) && use_precomputed
            @info "\t using precomputed solution gpe_3d.jld2" 
            gpe_3d = JLD2.load(join([save_path, join(["gpe_3d",add_string, ".jld2"])]))["gpe_3d"]
        else
            sol = runsim(sim_gpe_3d; info=true)
            gpe_3d = sol.u
            JLD2.save(join([save_path, join(["gpe_3d",add_string, ".jld2"])]), "gpe_3d",  gpe_3d)
        end
        # linear interpolation
        gpe_3d = sim_gpe_3d.psi_0
        x_axis = sim_npse.X[1] |> real
        x_axis_3d = sim_gpe_3d.X[1] |> real
        dx = sim_gpe_3d.X[1][2]-sim_gpe_3d.X[1][1]
        final_axial = Array(sum(abs2.(xspace(gpe_3d, sim_gpe_3d)), dims=(2, 3)))[:, 1, 1] * sim_gpe_3d.dV / dx |> real
        # we need to renormalize (error in the sum??)
        final_axial = final_axial / sum(final_axial * dx) |> real
        x_3d_range = range(-sim_gpe_3d.L[1]/2, sim_gpe_3d.L[1]/2, length(sim_gpe_3d.X[1])) 
        solution_3d = LinearInterpolation(x_3d_range, final_axial, extrapolation_bc = Line())
        plot!(p, x_axis, solution_3d(x_axis), label="GPE_3D", color=:red) 
        # q = plot(x_axis_3d, final_axial, label="GPE_3D", color=:red, linestyle=:dot) 
        # display(q)
        display(p)

        # s2 = estimate_sigma2(kspace(initial_3d, sim_gpe_3d), sim_gpe_3d)
        # sigma_2 = plot(x_axis_3d, s2, label="sigma2", color=:red, linestyle=:dot)
        # dens = sum(abs2.(initial_3d), dims=(2, 3))[:, 1, 1] * sim_gpe_3d.dV / dx |> real
        # plot!(sigma_2, x_axis_3d, dens, label="psi^2", color=:red)
        # display(sigma_2)


        # s2 = estimate_sigma2(gpe_3d, sim_gpe_3d
        # sigma_2 = plot(x_axis_3d, s2, label="sigma2", color=:red)
        # plot!(sigma_2, x_axis_3d, sigma2_old, label="NPSE", color=:red, linestyle=:dash)
        # plot!(sigma_2, x_axis_3d, sigma2_new, label="NPSE:plus", color=:red, linestyle=:dot)
        # plot!(sigma_2, x_axis_3d, final_axial, label="psi^2", color=:red)
        # display(sigma_2)

        # heatmap(abs2.(xspace(gpe_3d, sim_gpe_3d))[3, :, :], aspect_ratio=1, color=:viridis, title="GPE_3D")
        # savefig(p, "media/comparison_g$(gamma_param).png")

        # TODO implement and return also sigma2
        return p
    end
end