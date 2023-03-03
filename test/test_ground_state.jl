using CondensateDynamics 
using CUDA
## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (4.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
iswitch=-im
g = 0.0
V(x, t) = 1/2 * (x^2)
equation = GPE_1D
x = X[1]
@. psi_0 = 1/(pi^1/4) * exp(-x^2/20) # very broad initial state
psi_0 = psi_0/sqrt(norm_squared(psi_0, sim))
reltol=1e-7
@pack_Sim! sim

# Analytical solution: Gaussian
analytical_gs = zeros(N)
@. analytical_gs = exp(-(x^2)/2)/(pi^(1/4))

sol, err = testsim(sim)
@test err == false

numerical_gs = abs2.(sol[end])
@test isapprox(numerical_gs, abs2.(analytical_gs), rtol=1e-4)