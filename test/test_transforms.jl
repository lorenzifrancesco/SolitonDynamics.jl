using CondensateDynamics 
using CUDA
using OrdinaryDiffEq

L = (40.0,)
N = (1024,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)


@unpack_Sim sim
g = 0.0
equation = GPE_1D
iswitch = 1
x = X[1]
k = K[1]
dV= volume_element(L, N)
reltol = 1e-2

@. psi_0 = exp(-x^2/2/5) * exp(-im * x * 0.01)
psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))

psi_1 = exp.(-x.^2/2/5) .* exp.(-im * x * 0.01)
psi_1 = psi_1 / sqrt(sum(abs2.(psi_1) * dV))

initial_state = psi_0

alg = Tsit5()

@pack_Sim! sim

# inplace test for identity
kspace!(psi_0, sim)
xspace!(psi_0, sim)

@test isapprox(psi_0 - initial_state, zeros(N), atol=1e-9)


# out of place test for identity
psi_0 = initial_state
tmp = kspace(psi_0, sim)
tmp = xspace(tmp, sim)
x = LinRange(0, L[1], N[1]) |> collect

@test isapprox(tmp - psi_0, zeros(N), atol=1e-9)

# same result in place and out
psi_0 = initial_state
kspace!(psi_0, sim)
outofplace = xspace(psi_0, sim)
xspace!(psi_0, sim)

@test isapprox(outofplace - psi_0, zeros(N), atol=1e-9)

# inplace norm conservation

psi_0 = initial_state

kspace!(psi_0, sim)
xspace!(psi_0, sim)

@test isapprox(ns(psi_0, sim) - ns(initial_state, sim), 0.0, atol=1e-9)

# test norms (Parseval's theorem)

@info ns(psi_1, sim)
@info nsk(kspace(psi_1, sim), sim)

@test isapprox(ns(psi_1, sim), nsk(kspace(initial_state, sim), sim), atol=1e-9)

# test the correctness of transforms in a known case 


# (gaussian wavepacket)
psi_constant = ones(sim.N) |> complex
display(psi_constant)
psi_constant_trans = kspace(psi_constant, sim)
display(psi_constant_trans)