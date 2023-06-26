includet("CDResearch.jl")

function sdmf()
    if Threads.nthreads() == 1
        @warn "running in single thread mode!"
    else
        @info "running in multi-thread mode: n_threads =" Threads.nthreads()
    end

    ## procedure:
    # 1- compute ground states
    # 2- prepare dynamical simulations using GS
    # 3- compute tiles
    # 4- compute lines
    save_path = "results/"
    gamma_list = [0.15, 0.4, 0.6]
    for gamma in gamma_list
        @info "==== Using gamma: " gamma
        sd = load_parameters(gamma_param=gamma)
        @info "Required simulations: " keys(sd)
        # prepare ground states (saving them)
        if isfile(save_path * "gs_dict.jld2")
            @info "Loading GS library..."
            gs_dict = JLD2.load(save_path * "gs_dict.jld2")
        else
            @info "No GS library found! Saving an empty one..."
            gs_dict = Dict()
            JLD2.save(save_path * "gs_dict.jld2", gs_dict)
        end

        # preparing all simulations
        for (name, sim) in sd
            if haskey(gs_dict, hs(name, gamma))
                @info "Found in library item " (name, gamma)
            else
                @info "Computing item " (name, gamma)
                uu = get_ground_state(sim)
                push!(gs_dict, hs(name, gamma) => uu)
                JLD2.@save(save_path * "gs_dict.jld2", "sd", gs_dict)
            end
        end


        # visualize the ground states
        print("\nProceed to tiling? [Y/n]: ")
        if readline()[1] == 'n'
            @info "Aborting..."
            return
        end
        # run the tiling
        for (name, sim) in sd
            @info "Tiling " name
            tile = get_tiles(sim, name)
            JLD2.@save(save_path * "tile_dict.jld2", hs(name, gamma) => tile)
        end
    end
end