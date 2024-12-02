using Printf
function optical_lattice(v0, d, space)
  # See Eq.(3) of [PRA 75 033622 (2007)]
  # notice that the k_L is different from the one of Stratclyde
  return -v0 * cos.(2 * pi / d * space)
end

print("\n"*"@"^20*"BEGIN OF RUN"*"@"^20*"\n")
cmd = `date`
print(read(cmd, String))
using CSV, DataFrames, Tables, Printf, FFTW
using SolitonDynamics, Plots; gr()
print("\n@@@ packages are loaded ")
("\n"*"="^50*"n")

# params = [(0, 1), (1, 1)] # (p, S)
params = [(0, 0)]
# params = [(0, 1), (1, 1)]
# params = [(j, i) for i in range(0, 6) for j in range(1, 1)]
E = 1 # pulse energy or number of particles
gamma = 0.5
N = 200
L = 20.0
x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
sim = init_sim((L,), (N,))
sim.iswitch = 1.0
y_values = []
# params = range(-10, 10, 3)
print("\n p  S   gpS\n")
for (idx, par) in enumerate(params)
  sim.p = 0.0
  sim.S = 0.0
  sim.xi = (2*sim.p + sim.S + 1)
  # gamma_base_normalizzata = 4/3-0.05
  # gamma_base_normalizzata = 3.5 / 2
  sim.g = gamma2g(gamma, sim)
  sim.V0 = optical_lattice(1.8, 1, sim.X[1])
  sim.sigma2 = init_sigma2(sim.g)
  println("::::::::::", sim.g)
  @info "sim.g="*string(sim.g)
  sim.reltol = 1e-9
  sim.abstol = 1e-9
  sim.psi_0 = kspace(complex(gaussian(x/2, sim)), sim)
  sim.maxiters = 1e4
  sim.Nt = 100
  sim.tf = 0.5
  sim.t = LinRange(sim.ti, sim.tf, sim.Nt)
  sol = runsim(sim, info=true)
  print(sol)

  psi2 = Array{Float64,2}(undef, (length(sol[1].u), length(sol[1].u[1])))
  for (ix, u_t) in enumerate(sol[1].u)
    psi2[ix, :] = abs2.(xspace(u_t, sim))
  end
  # print(psi2)
  print(size(psi2))
  output = [vcat(sim.t[i], psi2[i, :]) for i in 1:size(psi2, 1)]
  CSV.write("results/experiment"*string(idx)*".csv", DataFrame(output, :auto))
  # sigma = sqrt.(sqrt.(1 .+ sim.g*psi2))
  # CSV.write("results/sigma"*string(idx)*".csv", DataFrame(x = x, y = sigma))
end

# Convert 2D array to a vector of vectors
# y_values1 = abs2.(gpe_analytical.(x, g2gamma(-6.0, NPSE)))
# CSV.write("results/axial3.csv", DataFrame(x = x, y = y_values1))

print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n"*"@"^20*"END OF RUN"*"@"^20*"\n")