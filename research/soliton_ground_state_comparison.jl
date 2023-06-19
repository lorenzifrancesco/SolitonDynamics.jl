```
 max g allowable for hashing = -5.0, 5.0
```
function hs(eq::String, g::Float64)
    @assert eq in ["G1", "N", "Np", "G3"]
    n = 0 
    if eq == "G1"
        n += 0
    elseif eq == "N"
        n += 1000
    elseif eq == "Np"
        n += 2000
    else
        n += 3000
    end
    n += Int(round(g * 100))
    @info n
    return n
end

function ihs(n::Int)
    if n < 500
        return ("G1", n/1000)
    elseif n < 1500
        return ("N", (n-1000)/1000)
    elseif n < 2500
        return ("Np", (n-2000)/1000)
    else
        return ("G3", (n-3000)/100)
    end
end

function all_ground_states()
    sd = load_parameters_gs()
    @assert all([s.iswitch for s in values(sd)] .== -im)
    save_path = "results/"

    if isfile("gs_dict.jld2")
        gs_dict = JLD2.load(join([save_path, "gs_dict.jld2"]))
        @info "Found GS dictionary: " gs_dict
    else
        gs_dict = Dict{Int64, AbstractArray}()
        # test save
        JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)
    end

    gamma_param_list = [0.15, 0.4, 0.6]
    use_precomputed = false
    take_advantage = true
    @info "Starting simulations..."
    for gamma_param in gamma_param_list
        # update simulation parameters
        @warn keys(sd)
        set_g!.(values(sd), gamma_param)
        sim_gpe_1d =    sd["G1"]
        sim_npse =      sd["N"]
        sim_npse_plus = sd["Np"]
        sim_gpe_3d =    sd["G3"]
        # =========================================================
        Plots.CURRENT_PLOT.nullableplot = nothing
        x = sim_gpe_1d.X[1] |> real
        analytical_gs = zeros(sim_gpe_1d.N[1])
        @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param))
        p = plot_final_density([analytical_gs], sim_gpe_1d; label="analytical", color=:orange, doifft=false, ls=:dashdot, title="gamma = $gamma_param")
        # plot_final_density!(p, [sim_], sim_gpe_1d; label="initial_GPE_1D", color=:grey, doifft=false, ls=:dashdot)

        if haskey(gs_dict, hs("G1", gamma_param))
            if use_precomputed
                @info "\t using precomputed solution G1" 
            else
                @info "\t deleting and recomputing solution G1"
                delete!(gs_dict, hs("G1", gamma_param))
                sol = runsim(sim_gpe_1d; info=true)
                push!(gs_dict, hs("G1", gamma_param) => sol.u)
            end
        else
            @info "computing G1" 
            sol = runsim(sim_gpe_1d; info=true)
            push!(gs_dict, hs("G1", gamma_param) => sol.u)
        end
        gpe_1d = gs_dict[hs("G1", gamma_param)]
        plot_final_density!(p, [gpe_1d], sim_gpe_1d; label="GPE_1D", color=:blue, ls=:dash)
        @warn keys(gs_dict)
        JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)

        # estimate width
        initial_sigma_improved = sqrt(sum(abs2.(sim_gpe_1d.X[1])  .* abs2.(gpe_1d) * sim_gpe_1d.dV) |> real)
        @warn "initial sigma improved" initial_sigma_improved
        if take_advantage 
            sim_npse.psi_0 = gpe_1d
        end
        if haskey(gs_dict, hs("N", gamma_param))
            if use_precomputed
                @info "\t using precomputed solution N" 
            else
                @info "\t deleting and recomputing solution N"
                delete!(gs_dict, hs("N", gamma_param))
                sol = runsim(sim_npse; info=true)
                push!(gs_dict, hs("N", gamma_param) => sol.u)
            end
        else
            @info "computing N"
            sol = runsim(sim_npse; info=true)
            push!(gs_dict, hs("N", gamma_param) => sol.u)
        end
        npse = gs_dict[hs("N", gamma_param)]
        plot_final_density!(p, [npse], sim_npse; label="NPSE", color=:green, ls=:dot)
        JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)

        if take_advantage
            sim_npse_plus.psi_0 = npse
        end
        if haskey(gs_dict, hs("N", gamma_param))
            if use_precomputed
                @info "\t using precomputed solution Np" 
            else
                @info "\t deleting and recomputing solution Np"
                delete!(gs_dict, hs("Np", gamma_param))
                sol = runsim(sim_npse_plus; info=true)
                push!(gs_dict, hs("Np", gamma_param) => sol.u)
            end
        else
            @info "computing Np" 
            sol = runsim(sim_npse_plus; info=true)
            push!(gs_dict, hs("Np", gamma_param) => sol.u)
        end
        npse_plus = gs_dict[hs("Np", gamma_param)]
        plot_final_density!(p, [npse_plus], sim_npse_plus; label="NPSE_der", ls=:dash, color=:green)
        JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)

        if take_advantage
            x = Array(sim_gpe_3d.X[1])
            y = Array(sim_gpe_3d.X[2])
            z = Array(sim_gpe_3d.X[3])
            initial_sigma_improved *= 5
            tmp = [exp(-((x/initial_sigma_improved)^2 + (y^2 + z^2)/2)) for x in x, y in y, z in z]
            sim_gpe_3d.psi_0 = CuArray(tmp)
            sim_gpe_3d.psi_0 .= sim_gpe_3d.psi_0 / sqrt(sum(abs2.(sim_gpe_3d.psi_0) * sim_gpe_3d.dV))
            initial_3d = copy(sim_gpe_3d.psi_0)
            kspace!(sim_gpe_3d.psi_0, sim_gpe_3d)
        end
        if haskey(gs_dict, hs("G3", gamma_param))
            if use_precomputed
                @info "\t using precomputed solution G3" 
            else
                @info "\t deleting and recomputing solution G3"
                delete!(gs_dict, hs("G3", gamma_param))
                sol = runsim(sim_gpe_3d; info=true)
                push!(gs_dict, hs("G3", gamma_param) => sol.u)
            end
        else
            @info "computing GPE_3D" 
            sol = runsim(sim_gpe_3d; info=true)
            push!(gs_dict, hs("G3", gamma_param) => sol.u)
        end
        gpe_3d = gs_dict[hs("G3", gamma_param)]
        JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)

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
            display(p)
        return p
    end
end