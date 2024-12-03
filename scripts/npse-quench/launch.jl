using Printf
using TOML

print("\n" * "@"^20 * "BEGIN OF RUN" * "@"^20 * "\n")
cmd = `date`
print(read(cmd, String))
using CSV, DataFrames, Tables, Printf, FFTW
using SolitonDynamics, Plots;
gr();
print("\n@@@ packages are loaded ")
("\n" * "="^50 * "n")

# constants: all units are in SI
hbar_nostro = 1.0546e-34
a_0 = 5.292e-11

function optical_lattice(v0, d, space)
  # See Eq.(3) of [PRA 75 033622 (2007)]
  # notice that the k_L is different from the one of Stratclyde
  return -v0 * cos.(2 * pi / d * space)
end

begin
  cf = TOML.parsefile("input/config.toml")
  # params = [(0, 1), (1, 1)] # (p, S)
  params = [(0, 0)]
  E = 1 # pulse energy or number of particles
  l_perp = sqrt(hbar_nostro/(cf["omega_perp"]*cf["m"]))
  @info l_perp
  g = 4 * pi * hbar_nostro^2 * cf["a_s"] * a_0 * cf["n_atoms"] / (cf["m"] * 2 * pi * l_perp^2)
  gamma = - cf["n_atoms"] * cf["a_s"] * a_0 / l_perp
  @info gamma
  g = gamma2g(gamma, sim)
  N = cf["n"]
  L = cf["l"]
  g5 = cf["l_3"]
  println(g5)
  x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
  dx = x[2] - x[1]
  sim = init_sim((L,), (N,))
  sim.iswitch = 1.0
  y_values = []
  # params = range(-10, 10, 3)
  print("\n p  S   gpS\n")
  for (idx, par) in enumerate(params)
    sim.p = 0.0
    sim.S = 0.0
    sim.xi = (2 * sim.p + sim.S + 1)
    # gamma_base_normalizzata = 4/3-0.05
    # gamma_base_normalizzata = 3.5 / 2
    sim.g = g
    sim.g5 = g5
    sim.V0 = optical_lattice(1.8, cf["d"]/l_perp, sim.X[1])
    sim.sigma2 = init_sigma2(sim.g)
    println("::::::::::", sim.g)
    @info "sim.g=" * string(sim.g)
    sim.reltol = 1e-9
    sim.abstol = 1e-9
    sim.psi_0 = kspace(complex(gaussian(x / 1.8, sim)), sim)
    sim.maxiters = 1e4
    sim.Nt = 100
    sim.tf = 1
    sim.t = LinRange(sim.ti, sim.tf, sim.Nt)
    sim.time_steps = Int64(ceil((sim.tf-sim.ti)/sim.dt))
    sol = runsim(sim, info=true)
    print(sol)

    psi2 = Array{Float64,2}(undef, (length(sol[1].u), length(sol[1].u[1])))
    # px_psi = Array{Float64,2}(undef, (length(sol[1].u), length(sol[1].u[1])))
    for (ix, u_t) in enumerate(sol[1].u)
      # psi2[ix, :] = vcat(real(im * diff(xspace(u_t, sim)) / dx), 0.0)
      psi2[ix, :] = abs2.(xspace(u_t, sim))
    end
    # print(psi2)
    print(size(psi2))
    output = [vcat(sim.t[i], psi2[i, :]) for i in 1:size(psi2, 1)]
    CSV.write("results/experiment" * string(idx) * ".csv", DataFrame(output, :auto))
    # sigma = sqrt.(sqrt.(1 .+ sim.g*psi2))
    # CSV.write("results/sigma"*string(idx)*".csv", DataFrame(x = x, y = sigma))
  end
end
# Convert 2D array to a vector of vectors
# y_values1 = abs2.(gpe_analytical.(x, g2gamma(-6.0, NPSE)))
# CSV.write("results/axial3.csv", DataFrame(x = x, y = y_values1))

print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n" * "@"^20 * "END OF RUN" * "@"^20 * "\n")