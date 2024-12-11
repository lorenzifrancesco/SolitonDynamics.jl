using Printf
using TOML
using PyFormattedStrings

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


begin
  # GS phase diagram:

  cf = TOML.parsefile("input/phase_diagram.toml")
  e_recoil = (pi * hbar / cf["d"])^2 / cf["m"]
  l_perp = sqrt(hbar_nostro / (cf_pre_quench["omega_perp"] * cf_pre_quench["m"]))
  e_perp = hbar_nostro * cf_pre_quench["omega_perp"]
  t_perp = cf_pre_quench["omega_perp"]^(-1)
  @info "l_perp" l_perp
  N = cf["n"]
  L = cf["l"]
  x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
  dx = x[2] - x[1]
  sim = init_sim((L,), (N,))

  # Phase diagram setup

  a_s_list = LinRange(cf["a_s_min"] * a_0,      cf["a_s_max"] * a_0,      cf["n_a_s"])
  v_0_list = LinRange(cf["v_0_min"] * e_recoil, cf["v_0_max"] * e_recoil, cf["n_v_0"])

  open("results/phase_diagram.txt", "w") do file
    write(file, f"===v_0 in [{v_0_min:.3e}, {v_0_max:.3e}], a_s in[{a_s_min:.3e}, {a_s_max:.3e}]===\n")

    for (iv, v_0) in enumerate(v_0_list)
      for (ia, a_s) in enumerate(a_s_list)
        write(file, f" {ia:4d}, {iv:4d} |")
        sim.dt = 1e-3
        sim.g5 = 0.0
        v_0_norm = v_0 / e_perp
        gamma = - a_s * cf["n_atoms"] / l_perp
        sim.g = gamma2g(gamma, sim)
        @info sim.g
        println(f"gamma = {gamma:.3e}")
        sim.V0 = optical_lattice(v_0_norm,
          cf["d"] / l_perp,
          0.0,
          0.0,
          sim.X[1])
        sim.sigma2 = init_sigma2(sim.g)
        sim.reltol = 1e-6
        sim.abstol = 1e-6
        sim.psi_0 = kspace(complex(gaussian(x / 2, sim)), sim)
        sim.iswitch = -im

        # GS finding
        sim.maxiters = 1e6
        sol = runsim(sim, info=true)
        try
          sim.psi_0 = sol.u
          psi2 = abs2.(xspace(sol.u, sim))
        catch
        end
        psi2 = zeros(cf["n"])
        for i in psi2
          write(file, f" {i:10.9f}")
        end
        write(file, f"\n")
        flush(file)
      end
    end
  end
end

print("@"^50)
cmd = `date`
print(read(cmd, String))
print("\n" * "@"^20 * "END OF RUN" * "@"^20 * "\n")