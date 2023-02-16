import CondensateDynamics: V, GPE_1D, GPE_3D, NPSE 
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
@. psi_0 = 1/(pi^1/4) * exp(-x^2/20)
@pack_Sim

# Analytical solution: Gaussian
analytical_gs = zeros(N)
@. analytical_gs = 1/(pi^1/4) * exp(-x^2/2)

sol, err = testsim(sim)
@test err == false
display(abs2.(sol[end]))

numerical_gs = abs.(xspace!(sol[end], sim))
display(numerical_gs)
@test isapprox(numerical_gs, analytical_gs, rtol=1e-4)