using Printf
print("\n"*"@"^20*"BEGIN OF RUN"*"@"^20*"\n")
cmd = `date`
print(read(cmd, String))
# print("NAME: ")
# @printf("[%30s]", readline())

using CSV, DataFrames, Tables, Printf, FFTW
using SolitonDynamics, Plots; gr()

print("\n@@@ packages are loaded ")

# TODO load configuration from input
("\n"*"="^50*"n")
params = [(0, 1), (1, 1), (2, 1)] # (p, S)
for (idx, par) in enumerate(params)
  N = 200
  L = 10.0
  sim = init_sim((L,), (N,))
  sim.g = -4
  sim.p = par[1]
  sim.S = par[2]
  sol = runsim(sim)
  x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
  sim.psi_0 = kspace(complex(exp.(x.^2)), sim)
  print(sol)
  CSV.write("results/data"*string(idx)*".csv", DataFrame(x = x, y = fftshift(abs2.(xspace(sol[1].u[1], sim)))))
end
print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n"*"@"^20*"END OF RUN"*"@"^20*"\n")