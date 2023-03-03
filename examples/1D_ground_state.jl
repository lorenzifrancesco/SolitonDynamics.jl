using Makie, LaTeXStrings, GLMakie
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
#gr(colorbar=false,size=(600,150),legend=false,grid=false)

function ground_state_1D()
    ## Solve the 1D harmonic oscillator
    # problem with 1D-GPE 
    L = (10.0,)
    N = (256,)
    sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)


    @unpack_Sim sim
    g = 0.0
    equation = GPE_1D
    iswitch = -im
    x = X[1]
    dV= volume_element(L, N)
    reltol = 1e-2
    @. psi_0 = exp(-x^2/2/5)
    psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    initial_state = psi_0
    @info norm_squared(psi_0, sim)
    alg = Tsit5()

    ## TODO : not getting the potential because of function
    @. V0= 1/2 * (x^2)

    @pack_Sim! sim

    # Analytical solution: Gaussian
    analytical_gs = zeros(N)
    @. analytical_gs = exp(-(x^2)/2)/(pi^(1/4))
    @info norm_squared(analytical_gs, sim)

    sol = runsim(sim; info=false)
    plot= lines(real.(x), abs.(V.(x, 0.0)),linestyle = :dash, show=true)
    lines!(real.(x), abs2.(sol[end]), show=true)
    lines!(real.(x), abs2.(analytical_gs), color=:red)
    lines!(real.(x), abs2.(initial_state), color=:orange, linestyle=:dot)
    @info "final distribution: " norm_squared(sol[end], sim)
    display(plot)

    # fig, ax = lines()
    # nframes = sim.Nt
    # framerate = 30
    # hue_iterator = range(0, 360, length=nframes)
    # i=1
    # record(fig, "evolution.gif", hue_iterator;
    #         framerate = framerate) do 
    #             lines!(real.(x), abs2.(sol[i]))
    #     i += 1
    # end
    return nothing
end