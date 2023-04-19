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

## Solve the 1D ground state 
# problem with 1D-GPE 
L = (40.0,)
N = (2048,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

# =========================================================
@unpack_Sim sim
iswitch = -im
equation = NPSE
manual = true
solver = SplitStep 
time_steps = 10000
g = -1.5
g_param = abs(g) / (2)
@info "gamma: " g_param
if equation==NPSE && g_param > 2/3
    @warn "we should expect NPSE collapse"
end
mu_analytical = (1 - g_param^2/2)
abstol = 1e-6
maxiters = 30000
dt = 0.022


x = X[1]
k = K[1]
dV= volume_element(L, N)

flags = FFTW.EXHAUSTIVE

width = 7
#@. V0 = 1/2 * (x^2/(width^2))
tf = Inf
# SPR condensate bright soliton t in units of omega_perp^-1
analytical_gs = zeros(N)
@. analytical_gs = sqrt(g_param/2) * 2/(exp(g_param*x) + exp(-x*g_param))
psi_0 .= exp.(-(x/10).^2)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))
initial_state .= psi_0

kspace!(psi_0, sim)

@pack_Sim! sim

# =========================================================
# Analytical solution: Gaussian
width = sqrt(width)
#@. analytical_gs = exp(-x^2/((width^2)*2)) /(pi^(1/4)*sqrt(width))

sol = runsim(sim; info=true)
final = sol.u

@info "chempot of analytical" chempot(analytical_gs, sim)
@info "final chempot" chempotk(final, sim)
@info "analytical calculation of chempot" mu_analytical

p = plot(real.(k), abs2.(kspace(initial_state, sim)), color=:blue, ls=:dot, lw=3, label="initial")

plot!(p, real.(k), abs2.(final), color=:red, label="final")
plot!(p, real.(k), abs2.(kspace(analytical_gs, sim)), ls=:dot, lw=2, color=:grey, label="analytical")
display(p)

xspace!(final, sim)

middle = Int(round(N[1]/2))
reduced_x = real.(x)[middle-10:middle+10]
p = plot(real.(x), abs2.(initial_state), color=:green, ls=:dash, lw=2, label="initial")
plot!(p, real.(x), abs2.(final), color=:red, label="final")
plot!(p, real.(x), abs2.(analytical_gs), ls=:dot, lw=2, color=:black, label="analytical")
#plot!(p, real.(x), abs.(V0), ls=:dash, lw=2, color=:green, label="potential")

display(p)