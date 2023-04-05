using LaTeXStrings, Plots
using CUDA
import GR

function plot_axial_heatmap(u, time_axis, sim::Sim{1, Array{ComplexF64}}; info=false)
    @unpack t, X = sim; x = X[1]
    u = reduce(hcat, u)
    u = mapslices(x->xspace(x, sim),u,dims=(1))
    ht = heatmap(real.(x), time_axis, abs2.(u)')
    display(ht)
    return ht
end

function plot_axial_heatmap(u, time_axis, sim::Sim{3, CuArray{ComplexF64}}, axis; info=false)
    @unpack t, X = sim; x = Array(X[axis])
    @assert axis == 3
    ax_list = (1, 2, 3)
    ax_list= filter(x->x!=axis, ax_list)

    [xspace!(x, sim) for x in u]
    # SPECIALIZE to axis = 3
    u_axial = [sum(abs2.(x), dims=ax_list)[1,1,:] for x in u]
    u_axial = Array(reduce(hcat, u_axial))
    ht = heatmap(real.(x), time_axis, u_axial')
    display(ht)
    return ht
end

function plot_final_density(u, psi_0, sim::Sim{1, Array{ComplexF64}}; info=false)
    @unpack t, X = sim; x = Array(X[axis])
    final = u[end]
    xspace!(final, sim)
    xspace!(psi_0, sim)
    info && @info "final norm" ns(final, sim)
    p = plot(real.(x), abs2.(psi_0), label="initial")
    plot!(p, real.(x), abs2.(final), label="final")
    display(p)
    return p
end

function plot_final_density(u, psi_0, sim::Sim{3, CuArray{ComplexF64}}, axis; info=false)
    @unpack t, X = sim; x = Array(X[axis])
    ax_list = (1, 2, 3)
    ax_list= filter(x->x!=axis, ax_list)
    final = u[end]
    final = xspace(final, sim)
    info && @info "final norm" ns(final, sim)
    final_axial = Array(sum(abs2.(final), dims=ax_list))[1,1,:]
    xspace!(psi_0, sim)
    psi_0_axial = Array(sum(abs2.(psi_0), dims=ax_list))[1,1,:]
    p = plot(real.(x), psi_0_axial, label="initial")
    plot!(p, real.(x), final_axial, label="final")
    display(p)
    return p
end