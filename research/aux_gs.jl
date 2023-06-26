includet("plotter_all.jl")
using ColorSchemes
plotly()
# no saves

function go()
    gamma = 0.6
    sd = load_parameters_alt(gamma_param=gamma)
    gg = sd["G3"]

    @unpack_Sim gg
    x = X[1] |> real
    @pack_Sim! gg
    true_min = chempot(gpe_analytical.(x, gamma; x0=gg.L[1] / 4), gg)
    @warn true_min
    p = plot(x, abs2.(gpe_analytical.(x, gamma; x0=gg.L[1] / 4)), label="analytical", color=:black)
    plot!(p, x, abs2.(xspace(gg.psi_0, gg)), label="initial", color=:grey)
    nn = 2
    mm = 5
    meas = []
    pal = palette([:red, :blue], nn)
    for i in 1:nn
        @unpack_Sim gg
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