function single_shot_dynamics()
    sd = load_parameters_dy(vv = 1., bb = 1., eqs = ["G1"], Nsaves = 5)
    for (name, sim) in sd
        plot_axial_heatmap(runsim(sim).u, runsim(sim).t, sim)
    end
end