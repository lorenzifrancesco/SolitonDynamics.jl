using CondensateDynamics

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (40.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim

iswitch=-im
g = 0.0
equation = GPE_1D
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
numerical_gs = xspace(sol[1], sim)

@test isapprox(ns((numerical_gs-analytical_gs), sim), 0.0, atol=1e-5)