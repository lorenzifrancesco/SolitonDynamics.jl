includet("plotter_all.jl")
using ColorSchemes
plotly()
# no saves

function explore_collapse()
    gamma = 0.6
    sd = load_parameters_alt(gamma_param=gamma)
    gg = sd["N"]

    @unpack_Sim gg
    x = X[1] |> real
    @pack_Sim! gg
    nn = 10
    mm = 5
    meas = []
    pal = palette([:red, :blue], nn)
    gamma_list = LinRange(0.0, 1.0, nn)
    #     gamma_list = LinRange(0.0, 0.6, nn)
    for i in 1:nn
        @unpack_Sim gg
        g = -2 * gamma_list[nn]
        dt = 0.01
        abstol = 1e-4
        reltol = 1e-4
        @pack_Sim! gg
        sol = runsim(gg; info=true)
        plot!(p, x, abs2.(xspace(sol.u, gg)), color=pal[i])
        push!(meas, chempotk(sol.u, gg))
    end
    @pack_Sim
    q = plot(1:nn, meas)
    plot!(q, 1:nn, ones(nn) * true_min, label="true min", color=:red)
    display(meas)
    display(p)
    display(q)
    return
end

function pinpoint_collapse()
    gamma = 0.6
    sd = load_parameters_alt(gamma_param=gamma)
    gg = sd["N"]

    @unpack_Sim gg
    x = X[1] |> real
    dt = 0.01
    abstol = 1e-4
    reltol = 1e-4
    @pack_Sim! gg

    iters = 25
    pal = palette([:red, :blue], iters)
    # gamma_list = LinRange(0.0, 0.6, nn)
    # we run a bisection
    gplus = 1
    gminus = 0
    for i in 1:iters
        gmid = (gplus + gminus) / 2
        gg.g = -2 * gmid
        try
            sol = runsim(gg; info=true)
            @info "going right"
            gminus = gmid
        catch err
            if isa(err, NpseCollapse)
                @warn "collapse at $gmid"
                @info "going left"
                gplus = gmid
            else
                throw(err)
            end
        end
    end
    @info "collapse point" (gplus + gminus) / 2
    return (gplus + gminus) / 2
end