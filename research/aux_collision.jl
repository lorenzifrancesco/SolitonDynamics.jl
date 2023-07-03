using ColorSchemes
# no saves

function collide()
    gamma = 0.65
    sd = load_parameters_alt(gamma_param=gamma; nosaves=true)
    prepare_for_collision!(sd, gamma)
    sim = sd["G3"]

    imprint_vel_set_bar!(sim; vv=1.0, bb=1.0)

    nn = 30
    meas = []
    pal = palette([:red, :blue], nn)
    
    mask_refl = map(xx -> xx > 0, sim.X[1] |> real)
    mask_tran = map(xx -> xx < 0, sim.X[1] |> real)
    
    if false
        @unpack_Sim gg
        dt = 0.01
        abstol = 1e-10
        reltol = 1e-10
        @pack_Sim! gg
        sol = runsim(gg; info=false)
        plot!(p, x, abs2.(xspace(sol.u, gg)), color=pal[1])
        @warn "chempot relative error" (chempotk(sol.u, gg) - true_min) / true_min
        push!(meas, chempotk(sol.u, gg))
    else
        dt_range = LinRange(0.001, 1, nn)
        for i in 1:nn
            @unpack_Sim sim
            dt = dt_range[i]
            time_steps  = Int(ceil((tf-ti)/dt))
            @warn "dt = " dt
            @warn "tf = " tf
            # @info "total_time " time_steps * dt # FIXME there is nothing we can do, but tf should not be present into the solver routine
            @pack_Sim! sim

            @time sol = runsim(sim; info=true)
            # plot_axial_heatmap(sol.u, sim.t, sim; show=true, title="dt=$dt")
            # heatmap(abs2.(sol.u))
            # display(ht)
            # readline()

            final = xspace(sol.u[end], sim)
            tran = ns(final, sim, mask_tran)
            refl = ns(final, sim, mask_refl)
            @info "T" tran
            @assert isapprox(tran+refl, 1.0, atol=1e-3)
            push!(meas, tran)
            GC.gc()
        end
        p=plot(dt_range, meas, color=pal)
        display(p)
    end
end