includet("plotter_all.jl")
using ColorSchemes
# no saves

function explore_collapse()
    gamma = 0.6
    sd = load_parameters_collapse(gamma_param=gamma)
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

function pinpoint_collapse(; dynamical::Bool=false)
    gamma = 0.6
    sd = load_parameters_collapse(gamma_param=gamma)
    gg = sd["N"]

    @unpack_Sim gg
    x = X[1] |> real
    dt = 0.001
    abstol = 1e-4
    if dynamical
        tf = 1000.0
        iswitch = 1
        t = LinRange(ti, tf, Nt)
    end
    reltol = 1e-4
    display(psi_0)
    @pack_Sim! gg

    iters = 15
    pal = palette([:red, :blue], iters)
    # gamma_list = LinRange(0.0, 0.6, nn)
    # we run a bisection
    decimals = 3
    
    gplus = 1
    gminus = 0.5
    diff = 1.0
    prev_gmid = 0.0
    cnt = 0
    # while diff > 10.0^(-decimals) && cnt < 30
    for i in 1:iters
        gmid = (gplus + gminus) / 2
        @info "trying" gmid
        gg.g = -2 * gmid
        gg.sigma2 = initi_sigma2(gg.g)
        try
            sol = runsim(gg; info=true)
            sol = nothing
            gminus = gmid
        catch err
            if isa(err, NpseCollapse)
                @warn "collapse at $gmid"
                gplus = gmid
            else
                throw(err)
            end
        end
        diff = abs(gmid - prev_gmid)
        prev_gmid = gmid
        cnt += 1
    end
    @info "collapse point" (gplus + gminus) / 2
    return (gplus + gminus) / 2
end

function pinpoint_collapse_3()
    gamma = 0.6
    sd = load_parameters_alt(gamma_param=gamma)
    gg = sd["G3"]

    @unpack_Sim gg
    x = X[1] |> real
    dt = 0.01
    abstol = 1e-4
    reltol = 1e-4
    @pack_Sim! gg

    iters = 50
    pal = palette([:red, :blue], iters)
    # gamma_list = LinRange(0.0, 0.6, nn)
    # we run a bisection
    gplus = 0.8
    gminus = 0.55
    @assert false
    for i in 1:iters
        gmid = (gplus + gminus) / 2
        @info "trying" gmid
        gg.g = -2 * gmid
        sol = runsim(gg; info=false)
        if maximum(abs2.(sol.u)) > 1/gg.dV
            @warn "collapse at $gmid"
            gplus = gmid
        else
            gminus = gmid
        end
    end
    @info "collapse point" (gplus + gminus) / 2
    return (gplus + gminus) / 2
end