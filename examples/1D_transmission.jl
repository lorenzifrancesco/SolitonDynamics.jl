using Pkg
Pkg.activate(".")

using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V
gr()
GR.usecolorscheme(1)

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (50.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
g = -0.587
gamma = abs(g) / 2
equation = GPE_1D
iswitch = 1
x = X[1]
k = K[1]
dV= volume_element(L, N)
reltol = 1e-3
x0 = L[1] / 4
vv = 0.5
tf = 2*x0/vv
maxiters = 20000
#@. psi_0 = exp(-x^2/2) * exp(-im*x*10)
@. psi_0 = sqrt(gamma/2) * 2/(exp(gamma*(x-x0)) + exp(-(x-x0)*gamma)) * exp(-im*x*vv)

psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()
@. V0 = 1*exp(-(x/0.699)^2)

@pack_Sim! sim


sol = runsim(sim; info=false)
@info "time steps" size(sol.t)
final = sol[end]
xspace!(final, sim)
xspace!(psi_0, sim)
@info "final distribution norm squared: " ns(final, sim)

p = plot(real.(x), abs2.(psi_0), label="initial")
plot!(p, real.(x), abs2.(final), label="final")
display(p)
xspace!(initial_state, sim)


map(x -> xspace!(x, sim), sol[:])
ht = heatmap(real.(x), sol.t, abs2.(sol)')
display(ht)
