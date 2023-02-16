using Makie, LaTeXStrings, GLMakie
using CondensateDynamics
#gr(colorbar=false,size=(600,150),legend=false,grid=false)

function ground_state_1D()
    ## Solve the 1D harmonic oscillator
    # problem with 1D-GPE 
    L = (10.0,)
    N = (256,)
    sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

    @unpack_Sim sim
    g = 0.0
    V(x, t) = 1/2 * (x^2)
    equation = GPE_1D
    iswitch = -im
    x = X[1]
    dV= volume_element(L, N)

    @. psi_0 = 1/(pi^1/4) * exp(-x^2/2)
    psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    @info norm_squared(psi_0, sim)
    @pack_Sim! sim
    display(sum(abs2.(sim.psi_0)*dV))

    # Analytical solution: Gaussian
    analytical_gs = zeros(N)
    @. analytical_gs = 1/(pi^1/4) * exp(-x^2/2)

    sol = runsim(sim)
    display(sum(abs2.(sol[end])) * 10/256)
    plot= lines(real.(x), abs.(V.(x, 0.0)),linestyle = :dash, show=true)
    lines!(real.(x), abs2.(sol[end]), show=true)
    # display(plot)

    return nothing
end