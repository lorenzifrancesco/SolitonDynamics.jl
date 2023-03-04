using Pkg
Pkg.activate(".")

using CondensateDynamics, OrdinaryDiffEq, LSODA
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
import GR

gr()
GR.usecolorscheme(1)

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (40.0,40.0,40.0)
N = (256,256,256)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
g = 0.0
equation = GPE_1D
iswitch = 1
x = X[1]
y = X[2]
z = X[3]
dV= volume_element(L, N)
reltol = 1e-4
tf = 10.0
@. psi_0 = exp(-(x^2+y^2+z^2)/2) * exp(-im*x*10)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
fft!(psi_0)
alg = Tsit5()
@. V0 = 1/2 * (x^2+y^2+z^2)
V(x,y,z,t) = x^2
@pack_Sim! sim


#sol = runsim(sim; info=false)
final = sol[end]
ifft!(final)
ifft!(psi_0)
@info "final distribution norm squared: " ns(final, sim)

p = plot(real.(x), abs2.(psi_0), label="initial")
plot!(p, real.(x), abs2.(final), label="final")
display(p)
ifft!(initial_state)


map(x -> ifft!(x), sol[:])
ht = heatmap(real.(x), sol.t, abs2.(sol)')
display(ht)