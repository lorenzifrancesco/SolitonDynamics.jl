plotly(size=(500, 500))
function single_shot_dynamics(sim::Sim{1, Array{Complex{Float64}}})
    u = runsim(sim).u
    plot_axial_heatmap(u, runsim(sim).t, sim)
    return u
end

function show_psi_0(sim::Sim{1, Array{Complex{Float64}}})
    @unpack N, L, X, psi_0, V0 = sim; x = X[1] |> real
    xpsi_0 = xspace(psi_0, sim)
    p = plot(x, abs2.(xpsi_0), label = "density")
    plot!(p, x, abs.(V0), ls=:dot, color=:grey)
    display(p)
end

function show_psi_0(sim::Sim{3, CuArray{Complex{Float64}}})
    @unpack dV, N, L, X, psi_0 = sim; x = X[1] |> real
    dx = x[2]-x[1]
    xpsi_0 = xspace(psi_0, sim)
    axial_density = sum(abs2.(xpsi_0), dims = (2,3))[:, 1, 1] * dV / dx
    p = plot(x, axial_density, label = "density")
    axial_potential = sum(abs2.(V0), dims = (2,3))[:, 1, 1]
    plot!(p, x, axial_potential, ls=:dot, color=:grey)
    display(p)
end

function explain_sim(sim::Sim{1, Array{Complex{Float64}}})
    @unpack_Sim sim
    if iswitch == 1 
        print("\n=== DYNAMICAL SIMULATION ===\n")
    else
        print("\n=== STATICAL  SIMULATION ===\n")
    end
    print("Equation type: $(equation)\n")
    if manual
        print("\t => Manual\n")
    else
        print("\t => Automatic\n")
    end
    print("\n_________Discretization info_________\n")
    print("\t Spatial domain: $(L[1])\n")
    print("\t Spatial points: $(N[1])\n")
    if iswitch == 1
        print("\t Time domain: $(ti) to $(tf)\n")
        print("\t (Time saves: $(Nt))\n")
    else
        print("\t Max iters: $(maxiters)\n")
    end
    print("\n_________Physical info_______________\n")
    print("\t g = $(g)\n")
    print("\t mu = $(mu)\n")
    print("\t Potential maximum = $(maximum(abs.(V0))), minimum = $(minimum(abs.(V0)))\n")
    @pack_Sim
    return nothing
end