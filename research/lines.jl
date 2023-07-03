function all_lines(; use_precomputed_lines=false)
    if Threads.nthreads() == 1
        @warn "running in single thread mode!"
    else
        @info "running in multi-thread mode: n_threads =" Threads.nthreads()
    end

    save_path = "results/"
    gamma_list = [0.55]

    for gamma in gamma_list
        @info "==== Using gamma: " gamma
    
        sd = load_parameters_alt(gamma_param=gamma; eqs=["G1"], nosaves=true)
        @info "Required simulations: " keys(sd)

        prepare_for_collision!(sd, gamma)

        if isfile(save_path * "line_dict.jld2")
            @info "Loading Lines library..."
            line_dict = JLD2.load(save_path * "line_dict.jld2")
        else
            @info "No lines library found! Saving an empty one..."
            line_dict = Dict()
            JLD2.save(save_path * "line_dict.jld2", line_dict)
        end

        for (name, sim) in sd
            @info "Tiling " name
            if haskey(line_dict, hs(name, gamma)) && use_precomputed_lines
                @info "Already found line for " name, gamma
            else
                line = get_lines(
                    sim, 
                    name; 
                    lines=3,
                    sweep="vel",
                    points=100)
                push!(line_dict, hs(name, gamma) => line)
                JLD2.save(save_path * "line_dict.jld2", line_dict)
            end
        end
    end
end

# TODO: extremely slow velocities are not even interesting
function get_lines(
    sim::Sim{1,Array{Complex{Float64}}},
    name::String="noname";
    lines=2,
    sweep="vel",
    points=100
    )


    saveto = "../media/lines_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    # asymmetric matrix: 
    @assert sweep in ["vel", "bar"]
    if sweep == "vel"
        vel_list = LinRange(0.1, max_vel, points)
        bar_list = LinRange(0.1, max_bar, lines) # FIXME find a better way to do this 0.1->1.0
        x_axis = vel_list
        y_axis = bar_list
    elseif sweep == "bar"
        vel_list = LinRange(0.1, max_vel, lines)
        bar_list = LinRange(0.1, max_bar, points)
        x_axis = bar_list
        y_axis = vel_list
    end
    tran = Array{Float64,2}(undef, (lines, points))
    refl = Array{Float64,2}(undef, (lines, points))

    @info "Filling sim grid..."
    sgrid = Array{Sim,2}(undef, (lines, points))
    archetype = sim
    sgrid[1, 1] = archetype

    @time begin
        for (iy, y) in enumerate(y_axis)
            for (ix, x) in enumerate(x_axis)
                if sweep == "vel"
                    sgrid[iy, ix] = imprint_vel_set_bar(archetype; vv=x, bb=y)
                elseif sweep == "bar"
                    sgrid[iy, ix] = imprint_vel_set_bar(archetype; vv=y, bb=x)
                end
            end
        end
    end
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, sgrid[1, 1].X[1] |> real)
    mask_tran = map(xx -> xx < 0, sgrid[1, 1].X[1] |> real)

    @info "Running lining..."
    if sweep == "vel"
        print("Bar\n|\n|\n|____Vel")
    elseif sweep == "bar"
        print("Vel\n|\n|\n|____Bar")
    end

    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(y_axis), enumerate(x_axis))
    full_time = @elapsed for ((iy, y), (ix, x)) in ProgressBar(iter)
        sim = sgrid[iy, ix]
        collapse_occured = false
        sol = nothing
        try
            avg_iteration_time += @elapsed sol = runsim(sim; info=false)
        catch err
            if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
                collapse_occured = true
            else
                throw(err)
            end
        end
        # catch maxiters hit and set the transmission to zero
        if sim.manual == false
            if sol.retcode != ReturnCode.Success
                @info "Run complete, computing transmission..."
                @info "Detected solver failure"
                tran[iy, ix] = 0.0
                refl[iy, ix] = 1.0
                @info "T = " tran[iy, ix]
            else
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[iy, ix] = ns(final, sim, mask_tran)
                refl[iy, ix] = ns(final, sim, mask_refl)
                @info "T = " tran[iy, ix]
            end
        else
            if !collapse_occured
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[iy, ix] = ns(final, sim, mask_tran)
                refl[iy, ix] = ns(final, sim, mask_refl)
            else
                @info "Run complete, detected collapse..."
                tran[iy, ix] = NaN
            end
            @info "T = " tran[iy, ix]
        end
        if ! isapprox(tran[iy, ix]+refl[iy, ix], 1.0, atol=1e-5)
            @warn "T+R != 1.0"
        end
    end
    @info "Tiling time            = " full_time
    @info "Total time in solver   = " avg_iteration_time
    @info "Average iteration time = " avg_iteration_time / lines^2

    JLD2.@save("tran_$(name).jld2", tran)
    JLD2.@save("refl_$(name).jld2", refl)
    return tran
end

"""
in the 3D case we do not have sufficient GPU mem, so we go serially
"""
function get_lines(
    archetype::Sim{3,CuArray{Complex{Float64}}},
    name::String="noname";
    lines=2,
    sweep="vel",
    points=100
    )
    saveto = "../media/lines_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    # asymmetric matrix: 
    @assert sweep in ["vel", "bar"]
    if sweep == "vel"
        vel_list = LinRange(0.1, max_vel, points)
        bar_list = LinRange(0.1, max_bar, lines) # FIXME find a better way to do this 0.1->1.0
        x_axis = vel_list
        y_axis = bar_list
    elseif sweep == "bar"
        vel_list = LinRange(0.1, max_vel, lines)
        bar_list = LinRange(0.1, max_bar, points)
        x_axis = bar_list
        y_axis = vel_list
    end
    tran = Array{Float64,2}(undef, (lines, points))
    refl = Array{Float64,2}(undef, (lines, points))


    @info "Proceeding serially from the archetype..."
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
    mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(y_axis), enumerate(x_axis))
    full_time = @elapsed for ((iy, y), (ix, x)) in ProgressBar(iter)
        sim = deepcopy(archetype)
        collapse_occured = false
        if sweep == "vel"
            imprint_vel_set_bar!(sim; vv=x, bb=y)
        elseif sweep == "bar"
            imprint_vel_set_bar!(sim; vv=y, bb=x)
        end
        sol = nothing
        try
            avg_iteration_time += @elapsed sol = runsim(sim; info=false)
        catch err
            if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
                collapse_occured = true
            else
                throw(err)
            end
        end
        # catch maxiters hit and set the transmission to zero
        if sim.manual == false
            if sol.retcode != ReturnCode.Success
                @info "Run complete, computing transmission..."
                @info "Detected solver failure"
                tran[iy, ix] = 0.0
                refl[iy, ix] = 1.0
                @info "T = " tran[iy, ix]
            else
                if !collapse_occured
                    final = sol.u[end]
                    @info "Run complete, computing transmission..."
                    xspace!(final, sim)
                    tran[iy, ix] = ns(final, sim, mask_tran)
                    refl[iy, ix] = ns(final, sim, mask_refl)
                else
                    @info "Run complete, detected collapse..."
                    tran[iy, ix] = NaN
                end
                @info "T = " tran[iy, ix]
            end
        else
            if !collapse_occured
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[iy, ix] = ns(final, sim, mask_tran)
                refl[iy, ix] = ns(final, sim, mask_refl)
            else
                @info "Run complete, detected collapse..."
                tran[iy, ix] = NaN
                refl[iy, ix] = NaN
            end
            @info "T = " tran[iy, ix]
        end
        if ! isapprox(tran[iy, ix]+refl[iy, ix], 1.0, atol=1e-5)
            @warn "T+R != 1.0"
        end
    end
    @info "Tiling time            = " full_time
    @info "Total time in solver   = " avg_iteration_time
    @info "Average iteration time = " avg_iteration_time / lines^2

    JLD2.@save("tran_$(name).jld2", tran)
    JLD2.@save("refl_$(name).jld2", refl)
    norm_bar = bar_list / max_bar
    norm_vel = vel_list / max_vel
    return tran
    return tran
end

function view_all_lines(; sweep="vel")
    line_file = "results/line_dict.jld2"
    @assert isfile(line_file)
    ld = load(line_file)
    for (k, v) in ld
        @info "found" ihs(k)
        @warn "check the size"
        if sweep == "vel"
            p = plot(xlabel="velocity", ylabel="barrier", title=ihs(k))
        elseif sweep == "bar"
            p = plot(xlabel="barrier", ylabel="velocity", title=ihs(k))
        end
        for iy in 1:size(v)[1]
            plot!(p, v[iy, :], label=string(iy))
        end
        savefig(p, "media/lines_" * string(ihs(k)) * ".pdf")
        #  display(p)
    end
end

# TODO compare equations