JULIA_CUDA_SOFT_MEMORY_LIMIT = "95%"

function all_tiles(; use_precomputed_tiles=false)
    if Threads.nthreads() == 1
        @warn "running in single thread mode!"
    else
        @info "running in multi-thread mode: n_threads =" Threads.nthreads()
    end

    save_path = "results/"
    gamma_list = [0.65, 0.55, 0.4, 0.3, 0.15]

    for gamma in gamma_list
        @info "==== Using gamma: " gamma
        sd = load_parameters_alt(gamma_param=gamma; nosaves=true)
        @info "Required simulations: " keys(sd)

        prepare_for_collision!(sd, gamma)
        if isfile(save_path * "tile_dict.jld2")
            @info "Loading Tiles library..."
            tile_dict = JLD2.load(save_path * "tile_dict.jld2")
        else
            @info "No Tiles library found! Saving an empty one..."
            tile_dict = Dict()
            JLD2.save(save_path * "tile_dict.jld2", tile_dict)
        end

        for (name, sim) in sd
            @info "Tiling " name
            if haskey(tile_dict, hs(name, gamma)) && use_precomputed_tiles
                @info "Already found tile for " name, gamma
            else
                tile = get_tiles(sim, name)
                push!(tile_dict, hs(name, gamma) => tile)
                JLD2.save(save_path * "tile_dict.jld2", tile_dict)
            end
        end
    end
end

function get_tiles(
    sim::Sim{1,Array{Complex{Float64}}},
    name::String="noname";
    tiles=100)

    saveto = "../media/tiles_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64,2}(undef, (tiles, tiles))
    refl = Array{Float64,2}(undef, (tiles, tiles))

    @info "Filling sim grid..."
    sgrid = Array{Sim,2}(undef, (tiles, tiles))
    archetype = sim
    sgrid[1, 1] = archetype
    @time begin
        for (vx, vv) in enumerate(vel_list)
            for (bx, bb) in enumerate(bar_list)
                sgrid[bx, vx] = imprint_vel_set_bar(archetype; vv=vv, bb=bb)
            end
        end
    end
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, sgrid[1, 1].X[1] |> real)
    mask_tran = map(xx -> xx < 0, sgrid[1, 1].X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
    full_time = @elapsed for ((vx, vv), (bx, bb)) in ProgressBar(iter)
        sim = sgrid[bx, vx]
        collapse_occured = false
        @info "Computing tile" (vv, bb)
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
                tran[bx, vx] = 0.0
                refl[bx, vx] = 1.0
                @info "T = " tran[bx, vx]
            else
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[bx, vx] = ns(final, sim, mask_tran)
                refl[bx, vx] = ns(final, sim, mask_refl)
                @info "T = " tran[bx, vx]
            end
        else
            if !collapse_occured
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[bx, vx] = ns(final, sim, mask_tran)
                refl[bx, vx] = ns(final, sim, mask_refl)
            else
                @info "Run complete, detected collapse..."
                tran[bx, vx] = NaN
            end
            @info "T = " tran[bx, vx]
        end
        if ! isapprox(tran[bx, vx]+refl[bx, vx], 1.0, atol=1e-5)
            @warn "T+R != 1.0"
        end
    end
    @info "Tiling time            = " full_time
    @info "Total time in solver   = " avg_iteration_time
    @info "Average iteration time = " avg_iteration_time / tiles^2

    JLD2.@save("tran_$(name).jld2", tran)
    JLD2.@save("refl_$(name).jld2", refl)
    norm_bar = bar_list / max_bar
    norm_vel = vel_list / max_vel
    return tran
end

"""
in the 3D case we do not have sufficient GPU mem, so we go serially
"""
function get_tiles(
    archetype::Sim{3,CuArray{Complex{Float64}}},
    name::String="noname";
    tiles=100)
    saveto = "../media/tiles_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64,2}(undef, (tiles, tiles))
    refl = Array{Float64,2}(undef, (tiles, tiles))

    @info "Proceeding serially from the archetype..."
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
    mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
    full_time = @elapsed for ((vx, vv), (bx, bb)) in ProgressBar(iter)
        sim = deepcopy(archetype)
        collapse_occured = false
        imprint_vel_set_bar!(sim; vv=vv, bb=bb)
        @info "Computing tile" (vv, bb)
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
                tran[bx, vx] = 0.0
                refl[bx, vx] = 1.0
                @info "T = " tran[bx, vx]
            else
                if !collapse_occured
                    final = sol.u[end]
                    @info "Run complete, computing transmission..."
                    xspace!(final, sim)
                    tran[bx, vx] = ns(final, sim, mask_tran)
                    refl[bx, vx] = ns(final, sim, mask_refl)
                else
                    @info "Run complete, detected collapse..."
                    tran[bx, vx] = NaN
                end
                @info "T = " tran[bx, vx]
            end
        else
            if !collapse_occured
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[bx, vx] = ns(final, sim, mask_tran)
                refl[bx, vx] = ns(final, sim, mask_refl)
            else
                @info "Run complete, detected collapse..."
                tran[bx, vx] = NaN
                refl[bx, vx] = NaN
            end
            @info "T = " tran[bx, vx]
        end
        if ! isapprox(tran[bx, vx]+refl[bx, vx], 1.0, atol=1e-5)
            @warn "T+R != 1.0"
        end
    end
    @info "Tiling time            = " full_time
    @info "Total time in solver   = " avg_iteration_time
    @info "Average iteration time = " avg_iteration_time / tiles^2

    JLD2.@save("tran_$(name).jld2", tran)
    JLD2.@save("refl_$(name).jld2", refl)
    norm_bar = bar_list / max_bar
    norm_vel = vel_list / max_vel
    return tran
    return tran
end





function view_all_tiles()
    tile_file = "results/tile_dict.jld2"
    @assert isfile(tile_file)
    td = load(tile_file)
    for (k, v) in td
        @info "found" ihs(k)
        axis = LinRange(0.0, 1.0, size(v)[1])
        ht = contour(axis, axis, v, clabels=true, xlabel="v", ylabel="b")
        savefig(ht, "media/tiles_" * string(ihs(k)) * ".pdf")
        #  display(ht)
    end
end