## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (40.0,)
N = (512,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim

iswitch=-im
g = 0.0
equation = GPE_1D
manual = true
abstol = 1e-7
x = X[1]
dV= volume_element(L, N)
psi_0 = exp.(-x.^2/5)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))

kspace!(psi_0, sim)

@. V0= 1/2 * (x^2)

@pack_Sim! sim

# Analytical solution: Gaussian
analytical_gs = zeros(N)
@. analytical_gs = exp(-(x^2)/2)/(pi^(1/4))

sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u, sim)
@test isapprox(ns((numerical_gs-analytical_gs), sim), 0.0, atol=1e-5)


## Solving the SPR-like soliton ground state
# problem with 1D-GPE 
L = (40.0,)
N = (4096,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N) 

@unpack_Sim sim

iswitch=-im
g_param = 0.6
g = - 2 * g_param

equation = GPE_1D
manual = true
abstol = 1e-12
x = X[1]
dV= volume_element(L, N)
psi_0 = exp.(-x.^2/5)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))

kspace!(psi_0, sim)

@pack_Sim! sim

# Analytical solution: soliton in 1D-GPE
analytical_gs = zeros(N)
@. analytical_gs = sqrt(g_param/2) * 2/(exp(g_param*x) + exp(-x*g_param))

sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u, sim)
@warn "too loose nonlinear 1D soliton test"
@info ns(numerical_gs, sim)
@info ns(analytical_gs, sim)
save("analytical_gs.jld2"; analytical_gs)
save("numerical_gs.jld2"; numerical_gs)
@test isapprox(ns((numerical_gs-analytical_gs), sim), 0.0, atol=1e-6)