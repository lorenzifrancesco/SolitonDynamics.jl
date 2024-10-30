using Printf

print("\n"*"@"^20*"BEGIN OF RUN"*"@"^20*"\n")
cmd = `date`
print(read(cmd, String))
using CSV, DataFrames, Tables, Printf, FFTW
using SolitonDynamics, Plots; gr()
print("\n@@@ packages are loaded ")
("\n"*"="^50*"n")


params = [(0, 1), (1, 1), (2, 1)] # (p, S)
params = [(0, 1), (1, 1)]
E = 1 # pulse energy or number of particles
N = 200
L = 10.0
x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
sim = init_sim((L,), (N,))
y_values = []
params = range(-10, 10, 3)
for (idx, par) in enumerate(params)
  sim.g = -(2 * (2/3) + par/10) * E # in SolitonDynamics the field is normalized to 1
  # sim.p = par[1]
  # sim.S = par[2]
  sim.reltol = 1e-9
  sim.abstol = 1e-9
  sim.psi_0 = kspace(complex(gaussian(x, sim)), sim)
  sim.maxiters = 1e6
  sol = runsim(sim, info=true)
  y_values = abs2.(xspace(sol.u[1], sim))
  CSV.write("results/data"*string(idx % 3 + 1)*".csv", DataFrame(x = x, y = y_values))
end

y_values1 = abs2.(gpe_analytical.(x, g2gamma(-6.0, NPSE)))
CSV.write("results/data3.csv", DataFrame(x = x, y = y_values1))

print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n"*"@"^20*"END OF RUN"*"@"^20*"\n")