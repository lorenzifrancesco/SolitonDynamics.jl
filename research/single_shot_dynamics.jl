sd = load_parameters_dy()
for (name, sim) in sd
    plot_axial_heatmap(runsim(sim).u, runsim(sim).t, sim)
end