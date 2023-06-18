JULIA_CUDA_SOFT_MEMORY_LIMIT ="95%"

# XXX remark: good idea to vectorize on equations
function get_tiles(eq = "G1")
    saveto = "../media/tiles_$(eq).pdf"
    @info "Setting tiles configuration..."    
    tiles = 10
    barrier_width = 0.5 
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64, 2}(undef, (tiles, tiles))
    refl = Array{Float64, 2}(undef, (tiles, tiles))

    @info "Loading parameters, filling sim grid..."
    sgrid = Array{Sim, 2}(undef, (tiles, tiles))
    archetype = prepare_in_ground_state(load_parameters_dy(eqs=[eq], Nsaves=2)[eq])
    sgrid[1, 1] = archetype
    @time begin
        for (vx, vv) in enumerate(vel_list)
            for (bx, bb) in enumerate(bar_list)
                sgrid[bx, vx] = imprint_vel_set_bar(archetype, vv, bb)
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

    JLD2.@save("tran.jld2", tran)
    JLD2.@save("refl.jld2", refl)
    norm_bar = bar_list / max_bar
    norm_vel = vel_list / max_vel
    ht = Plots.heatmap(norm_bar, norm_vel, tran')
    display(ht)
    Plots.savefig(ht, saveto)

    # display(sol.destats)
    # display(sol.retcode)
end

function all_tiles()
    return
end