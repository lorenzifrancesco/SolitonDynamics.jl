plotly(size=(800, 600))
function single_shot_dynamics(sim::Sim{1, Array{Complex{Float64}}})
    plot_axial_heatmap(runsim(sim).u, runsim(sim).t, sim)
end

function show_psi_0(sim::Sim{1, Array{Complex{Float64}}})
    @unpack N, L, X, psi_0 = sim; x = X[1] |> real
    xpsi_0 = xspace(psi_0, sim)
    p = plot(x, abs2.(xpsi_0), label = "density")
    display(p)
end

function show_psi_0(sim::Sim{3, CuArray{Complex{Float64}}})
    @unpack dV, N, L, X, psi_0 = sim; x = X[1] |> real
    dx = x[2]-x[1]
    xpsi_0 = xspace(psi_0, sim)
    axial_density = sum(abs2.(xpsi_0), dims = (2,3))[:, 1, 1] * dV / dx
    p = plot(x, axial_density, label = "density")
    display(p)
end
