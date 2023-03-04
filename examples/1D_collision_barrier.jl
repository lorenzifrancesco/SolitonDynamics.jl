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
L = (40.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
g = -3.0
equation = GPE_1D
iswitch = 1
x = X[1]
k = K[1]
dV= volume_element(L, N)
reltol = 1e-4
tf = 10.0
x0 = -10
@. psi_0 = exp(-(x-x0)^2/2) * exp(im*x*2)

psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = Tsit5()
@. V0 = 4* exp(-10*x^2)

@pack_Sim! sim


sol = runsim(sim; info=false)
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