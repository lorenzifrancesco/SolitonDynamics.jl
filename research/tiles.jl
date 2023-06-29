JULIA_CUDA_SOFT_MEMORY_LIMIT ="95%"

# XXX remark: good idea to vectorize on equations
function get_tiles(sim::Sim{1, Array{Complex{Float64}}}, name::String="noname"; tiles=4)
    saveto = "../media/tiles_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64, 2}(undef, (tiles, tiles))
    refl = Array{Float64, 2}(undef, (tiles, tiles))

    @info "Filling sim grid..."
    sgrid = Array{Sim, 2}(undef, (tiles, tiles))
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
    mask_refl = map(xx -> xx>0, sgrid[1, 1].X[1] |> real)
    mask_tran = map(xx -> xx<0, sgrid[1, 1].X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
    full_time = @elapsed for ((vx, vv), (bx, bb)) in ProgressBar(iter)
        sim = sgrid[bx, vx]
        @info "Computing tile" (vv, bb)
        avg_iteration_time += @elapsed sol = runsim(sim; info=false)
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
            final = sol.u[end]
            @info "Run complete, computing transmission..."
            xspace!(final, sim)
            tran[bx, vx] = ns(final, sim, mask_tran)
            refl[bx, vx] = ns(final, sim, mask_refl)
            @info "T = " tran[bx, vx]
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

function get_tiles(archetype::Sim{3, CuArray{Complex{Float64}}}, name::String="noname"; tiles=4)
    saveto = "../media/tiles_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64, 2}(undef, (tiles, tiles))
    refl = Array{Float64, 2}(undef, (tiles, tiles))

    @info "Proceeding serially from the archetype..."
    # all sims have the same x
    mask_refl = map(xx -> xx>0, archetype.X[1] |> real)
    mask_tran = map(xx -> xx<0, archetype.X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
    full_time = @elapsed for ((vx, vv), (bx, bb)) in ProgressBar(iter)
        sim = deepcopy(archetype)
        collapse_occured = false
        imprint_vel_set_bar!(sim; vv=vv, bb=bb)
        @info "Computing tile" (vv, bb)
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
            end
                @info "T = " tran[bx, vx]
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

function all_tiles()
    @info "Load parameters... "
    sd = load_parameters_dy(Nsaves=2)
    for (name, sim) in sd
        @info ">>> Computing for $name:"
        @info "Preparing archetype in ground state..."
        prepare_in_ground_state!(sim)
        @info "Computing tiles..."
        get_tiles(sim, name)
    end
    return nothing
end