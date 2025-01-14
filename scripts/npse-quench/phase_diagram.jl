using Printf
using TOML
using PyFormattedStrings
using Base.Threads
using DelimitedFiles

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
e_r = 2.1798723611030e-18

function optical_lattice(v_0, d, tilt, l_x, space)
  # See Eq.(3) of [PRA 75 033622 (2007)]
  # v_0 is in normalized units (hbar omega_perp)
  # notice that the k_L is different from the one of Stratclyde
  if l_x == 0.0
    return -v_0 * cos.(2 * pi / d * space) + tilt * space
  else
    return -v_0 * cos.(2 * pi / d * space) + tilt * space + 1 / 2 * (space / l_x) .^ 2
  end
end

# data = readdlm("results/phase_diagram.txt", '|')
# indices = []
# data = data[2:end]
# for d in data
#   try
#     push!(indices, map(x -> parse(Int, x), split(d, ",")))
#   catch 
#   end
# end
# println(indices)
indices = []
begin
  # GS phase diagram:
  cf = TOML.parsefile("input/phase_diagram.toml")
  e_recoil = (pi * hbar / cf["d"])^2 / (2 * cf["m"])
  l_perp = sqrt(hbar_nostro / (cf["omega_perp"] * cf["m"]))
  e_perp = hbar_nostro * cf["omega_perp"]
  t_perp = cf["omega_perp"]^(-1)
  @info "l_perp" l_perp
  N = cf["n"]
  L = cf["l"]
  x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
  dx = x[2] - x[1]
  sim = init_sim((L,), (N,))

  sim.reltol = 1e-6
  sim.abstol = 1e-6
  sim.maxiters = 1e6
  sim.g5 = 0.0
  sim.dt = 1e-3
  # Phase diagram setup

  a_s_list = LinRange(cf["a_s_min"] * a_0,      cf["a_s_max"] * a_0,      cf["n_a_s"])
  v_0_list = LinRange(cf["v_0_min"] * e_recoil, cf["v_0_max"] * e_recoil, cf["n_v_0"])
  v_0_min = cf["v_0_min"]
  v_0_max = cf["v_0_max"]
  a_s_min = cf["a_s_min"]
  a_s_max = cf["a_s_max"]
  print(a_s_list)
  print(v_0_list)
  open("results/phase_diagram.txt", "a") do file
    first_flag = true
    for (iv, v_0) in enumerate(v_0_list)
      for (ia, a_s) in enumerate(a_s_list)
        if !([ia, iv] in indices) 
          write(file, f" {ia:4d}, {iv:4d} |")
          v_0_norm = v_0 / e_perp
          gamma = - a_s * cf["n_atoms"] / (l_perp)
          sim.g = gamma2g(gamma / 2, sim)
          sim.sigma2 = init_sigma2(sim.g)
          println(f">{ia:4d}, {iv:4d} | gamma = {gamma:.3e}, g = {sim.g:.3e}")
          sim.V0 = optical_lattice(v_0_norm,
            cf["d"] / l_perp,
            0.0,
            0.0,
            sim.X[1])
          if first_flag == true
            sim.psi_0 = kspace(complex(gaussian(x / 3, sim)), sim)
            first_flag = false
          end
          sim.iswitch = -im

          # GS finding
          sol = runsim(sim, info=true)
          psi2 = zeros(cf["n"])
          try
            sim.psi_0 = sol.u

            psi2 = abs2.(xspace(sol.u, sim))
          catch
            print("miss! printing zeros")
            sim.psi_0 = kspace(complex(gaussian(x / 3, sim)), sim)
          end
          for i in psi2
            write(file, f" {i:10.9f}")
          end
          write(file, f"\n")
          flush(file)
        else
          print("Skipping...")
        end
      end
    end
  end
end

print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n" * "@"^20 * "END OF RUN" * "@"^20 * "\n")