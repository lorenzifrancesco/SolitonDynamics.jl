## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
# using Plots, Test, SolitonDynamics#dev
# ^ used in REPL refining

L = (40.0,)
N = (512,)
sim = Sim{length(L),Array{Complex{Float64}}}(L = L, N = N)

@unpack_Sim sim
iswitch = -im
g = 0.0
equation = GPE_1D
manual = true
abstol = 1e-7
x = X[1]
dV = volume_element(L, N)
psi_0 = exp.(-x .^ 2 / 5)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))
kspace!(psi_0, sim)
@. V0 = 1 / 2 * (x^2)
@pack_Sim! sim

# Analytical solution: Gaussian
analytical_gs = zeros(N)
@. analytical_gs = exp(-(x^2) / 2) / (pi^(1 / 4))

sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u, sim)
q = plot(real(x), abs2.(numerical_gs), color=:blue, label="num")
plot!(q, real(x), abs2.(analytical_gs), color=:red, label="ana")
savefig(q, "media/checks/test_gauss.pdf")
@test isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol = 1e-5)


## Solving the SPR-like soliton ground state
# problem with 1D-GPE 
L = (40.0,)
N = (512,)
sim = Sim{length(L),Array{Complex{Float64}}}(L = L, N = N)

g_param = 0.65

@unpack_Sim sim
iswitch = -im
g = - 2* g_param
equation = GPE_1D
manual = true
abstol = 1e-12
maxiters = 100000 # in the research is set to 1e10
x = X[1]
dV = volume_element(L, N)
psi_0 = exp.(-x .^ 2 / 10)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))
kspace!(psi_0, sim)
@pack_Sim! sim

# Analytical solution: soliton in 1D-GPE
analytical_gs = zeros(N)
@. analytical_gs = sqrt(g_param / 2) * 2 / (exp(g_param * x) + exp(-x * g_param))

sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u, sim)
@warn "too loose nonlinear 1D soliton test"
@info ns(numerical_gs, sim)
@info ns(analytical_gs, sim)
# save("analytical_gs.jld2"; analytical_gs)
# save("numerical_gs.jld2"; numerical_gs)
p = plot(real(x), abs2.(numerical_gs), color=:blue, label="num")
plot!(p, real(x), abs2.(analytical_gs), color=:red, label="ana")
savefig(p, "media/checks/test.pdf")
@test isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol = 1e-6)
