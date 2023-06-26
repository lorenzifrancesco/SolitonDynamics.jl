includet("plotter_all.jl")
using ColorSchemes
plotly()
# no saves

function go()
    gamma = 0.6
    sd = load_parameters_alt(gamma_param=gamma)
    gg = sd["G1"]

    @unpack_Sim gg
    x = X[1] |> real
    @pack_Sim! gg
    true_min = chempot(gpe_analytical.(x, gamma; x0=gg.L[1] / 4), gg)
    @warn true_min
    # p = plot(x, abs2.(gpe_analytical.(x, gamma; x0=gg.L[1] / 4)), label="analytical", color=:black)
    # plot!(p, x, abs2.(xspace(gg.psi_0, gg)), label="initial", color=:grey)
    nn = 15
    mm = 5
    meas = []
    pal = palette([:red, :blue], nn)
    for i in 1:nn
        @unpack_Sim gg
        dt = 0.1 / (i * 10) # FIXME tutto indica ad un problema di conservazione fisica
        abstol = 1e-3
        reltol = 1e-3
        @pack_Sim! gg
        sol = runsim(gg; info=true)
        # plot!(p, x, abs2.(xspace(sol.u, gg)), color=pal[i])
        push!(meas, chempotk(sol.u, gg))
    end
    @pack_Sim
    p = plot(1:nn, meas)
    plot!(p, 1:nn, ones(nn) * true_min, label="true min", color=:red)
    display(meas)
    display(p)
    return
end