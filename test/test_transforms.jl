using CondensateDynamics 
using CUDA
using OrdinaryDiffEq

L = (10.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)


@unpack_Sim sim
g = 0.0
equation = GPE_1D
iswitch = 1
x = X[1]
dV= volume_element(L, N)
reltol = 1e-2

@. psi_0 = exp(-x^2/2/5) * exp(-im * x * 0.01)
psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))

initial_state = psi_0

alg = Tsit5()

@pack_Sim! sim

# inplace test
original = psi_0
kspace!(psi_0, sim)
xspace!(psi_0, sim)

@test isapprox(psi_0 - original, zeros(N), atol=1e-10)

# out of place test
tmp = kspace(psi_0, sim)
result = xspace(tmp, sim)

@test isapprox(result - psi_0, zeros(N), atol=1e-10)