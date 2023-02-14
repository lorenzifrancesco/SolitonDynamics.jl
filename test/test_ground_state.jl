import CondensateDynamics: V, GPE_1D, GPE_3D, NPSE 
using CUDA
## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (4.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim

g = 0.0
V(x, t) = 1/2 * (x^2)
equation = GPE_1D
x = X[1]
@pack_Sim

# Analytical solution: Gaussian

@. analytical_gs = 1/(pi^1/4) * exp(-x^2/2)

sol, err = testsim(sim)
@test err==false

print(sol)
numerical_gs = xspace(sol)

@test isapprox(numerical_gs, analytical_gs, rtol=1e-4)