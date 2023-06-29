includet("plotter_all.jl")
using ColorSchemes
# no saves


function go()
    gamma = 0.6
    sd = load_parameters_alt(gamma_param=gamma)
    gg = sd["G1"]

    @unpack_Sim gg
    x = X[1] |> real
    solver = SplitStep
    psi_0 .= gpe_analytical.(x, gamma; x0=gg.L[1] / 4)
    kspace!(psi_0, gg)
    @pack_Sim! gg

    true_min = chempot(gpe_analytical.(x, gamma; x0=gg.L[1] / 4), gg)
    @info true_min
    p = plot(x, abs2.(gpe_analytical.(x, gamma; x0=gg.L[1] / 4)), label="analytical", color=:black)
    plot!(p, x, abs2.(xspace(gg.psi_0, gg)), label="initial", color=:grey)
    nn = 2
    mm = 5
    meas = []
    pal = palette([:red, :blue], nn)

    if true
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
        for i in 1:nn
            @unpack_Sim gg
            dt = 0.01
            abstol = 1e-4
            reltol = 1e-4
            @pack_Sim! gg
            sol = runsim(gg; info=true)
            plot!(p, x, abs2.(xspace(sol.u, gg)), color=pal[i])
            push!(meas, chempotk(sol.u, gg))
            display(q)
        end
    end
    # display(meas)
    # display(p)
    return
end