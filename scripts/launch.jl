using Printf

print("\n"*"@"^20*"BEGIN OF RUN"*"@"^20*"\n")
cmd = `date`
print(read(cmd, String))
using CSV, DataFrames, Tables, Printf, FFTW
using SolitonDynamics, Plots; gr()
print("\n@@@ packages are loaded ")
("\n"*"="^50*"n")

params = [(0, 1), (1, 1)] # (p, S)
# params = [(0, 0)]
# params = [(0, 1), (1, 1)]
# params = [(j, i) for i in range(0, 6) for j in range(1, 1)]
E = 1 # pulse energy or number of particles
N = 200
L = 20.0
x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
sim = init_sim((L,), (N,))
y_values = []
# params = range(-10, 10, 3)
print("\n p  S   gpS\n")
for (idx, par) in enumerate(params)
  sim.p = par[1]
  sim.S = par[2]
  sim.xi = (2*sim.p + sim.S + 1)
  gamma_base_normalizzata = 4/3-0.05
  gamma_base_normalizzata = 3.5 / 2
  sim.g = gpS(-gamma_base_normalizzata*2, sim) # in SolitonDynamics the field is normalized to 1
  sim.sigma2 = init_sigma2(sim.g)
  println("::::::::::", sim.g)
  @info "sim.g="*string(sim.g)
  # print("\n")
  sim.reltol = 1e-9
  sim.abstol = 1e-9
  sim.psi_0 = kspace(complex(gaussian(x, sim)), sim)
  sim.maxiters = 1e4
  sol = runsim(sim, info=true)
  psi2 = abs2.(xspace(sol.u, sim))
  CSV.write("results/axial"*string(idx)*".csv", DataFrame(x = x, y = psi2))
  sigma = sqrt.(sqrt.(1 .+ sim.g*psi2))
  CSV.write("results/sigma"*string(idx)*".csv", DataFrame(x = x, y = sigma))
end

y_values1 = abs2.(gpe_analytical.(x, g2gamma(-6.0, NPSE)))
CSV.write("results/axial3.csv", DataFrame(x = x, y = y_values1))

print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n"*"@"^20*"END OF RUN"*"@"^20*"\n")