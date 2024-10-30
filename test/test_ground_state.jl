## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
# using Plots, Test, SolitonDynamics#dev
# ^ used in REPL refining

L = (40.0,)
N = (512,)
sim = Sim{length(L),Array{ComplexF64}}(L=L, N=N)

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
analytical_gs = gaussian(x, sim)

sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u[end], sim)
@test isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol=1e-5)


## Solving the SPR-like soliton ground state
# problem with 1D-GPE 
L = (40.0,)
N = (256,)
sim = Sim{length(L),Array{ComplexF64}}(L=L, N=N)
g_param = 0.65
@unpack_Sim sim
iswitch = -im
g = -2 * g_param
equation = GPE_1D
manual = true
abstol = 1e-3
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
@info "=== Computing 1D-GPE"
@time sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u[end], sim)
@test isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol=1e-6)

## We need to have loose resembling also with NPSE and  NPSE+
@unpack_Sim sim
equation = NPSE
@pack_Sim! sim
@info "=== Computing NPSE"
@time sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u[end], sim)
@test isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol=1e-2)
@test !isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol=1e-9)

@unpack_Sim sim
equation = NPSE_plus
@pack_Sim! sim
@info "=== Computing NPSE+"
@time sol, err = testsim(sim)
@test err == false
numerical_gs = xspace(sol.u[end], sim)
@test isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol=5e-2)
@test !isapprox(ns((numerical_gs - analytical_gs), sim), 0.0, atol=1e-9)
@warn "We set 5% tolerance on L2 norm for NPSE+ vs 1D-GPE analytical solution."