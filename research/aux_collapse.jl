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
    p = plot(x, abs2.(gpe_analytical.(x, gamma; x0 = gg.L[1]/4)), label="analytical", color=:black)
    plot!(p, x, abs2.(xspace(gg.psi_0, gg)), label="initial", color=:grey)
    nn = 3
    pal = palette([:red, :blue], nn)
    gamma_range = LinRange(0, 1, nn)
    for i in 1:nn
        @unpack_Sim gg
        iswitch = 1
        g = -2 * gamma_range[i]
        dt = 0.00001 
        abstol = 1e-4
        reltol = 1e-4
        @pack_Sim! gg
        sol = runsim(gg; info=true)
        plot!(p, x, abs2.(xspace(sol.u, gg)), color=pal[i])
    end
    @pack_Sim
    display(p)
    return
end