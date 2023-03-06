using Pkg
Pkg.activate(".")

using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
using LSODA
import CondensateDynamics.V
import FFTW
gr()
GR.usecolorscheme(1)

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (40.0,)
N = (512,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim
iswitch = -im
g = -0.0
equation = GPE_1D
x = X[1]
k = K[1]
dV= volume_element(L, N)
reltol = 1e-4
psi_0 = exp.(-x.^2/5)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))
initial_state .= psi_0
@info "norm squared initial state" ns(psi_0, sim)
flags = FFTW.EXHAUSTIVE
kspace!(psi_0, sim)
## TODO : not getting the potential because of function
width = 1
@. V0= 1/2 * ((x/width)^2)
tf = Inf
@pack_Sim! sim

# Analytical solution: Gaussian
analytical_gs = zeros(N)
@. analytical_gs = exp(-((x/width)^2)/2)/(pi^(1/4)*sqrt(width))
@warn ns(analytical_gs, sim)
sol = runsim(sim; info=false)

final = sol[end] |> collect
@info "chemical potential" chempot(final, sim)
p = plot(real.(k), abs2.(kspace(initial_state, sim)), color=:blue, ls=:dot, lw=3, label="initial")
plot!(p, real.(k), abs2.(final), color=:red, label="final")
plot!(p, real.(k), abs2.(kspace(analytical_gs, sim)), ls=:dot, lw=2, color=:grey, label="analytical")
display(p)

# vettor = ones(N[1])
# vettor = xspace(vettor, sim)
# display(vettor)
xspace!(final, sim)
@info "norm final " ns(final, sim)

middle = Int(round(N[1]/2))
reduced_x = real.(x)[middle-10:middle+10]
p = plot(real.(x), abs2.(initial_state), color=:blue, ls=:dot, lw=3, label="initial")
plot!(p, real.(x), abs2.(final), color=:red, label="final")
plot!(p, real.(x), abs2.(final)/ns(final, sim), color=:red,  label="final normalized", ls=:dash)
plot!(p, real.(x), abs2.(analytical_gs), ls=:dot, lw=2, color=:grey, label="analytical")
display(p)

@warn sum(abs2.(final-analytical_gs))