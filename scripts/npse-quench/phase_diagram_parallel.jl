using Printf
using TOML
using PyFormattedStrings
using Base.Threads 

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

# Create a channel to collect results safel
function parallel_run()
    # GS phase diagram:
    cf = TOML.parsefile("input/phase_diagram.toml")
    e_recoil = (pi * hbar / cf["d"])^2 / cf["m"]
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
    a_s_list = LinRange(cf["a_s_min"] * a_0, cf["a_s_max"] * a_0, cf["n_a_s"])
    v_0_list = LinRange(cf["v_0_min"] * e_recoil, cf["v_0_max"] * e_recoil, cf["n_v_0"])

    # Open file and write the header
    open("results/phase_diagram.txt", "w") do file
        write(file, "===v_0 in [$(cf["v_0_min"]):.3e}, $(cf["v_0_max"]):.3e}], a_s in[$(cf["a_s_min"]):.3e}, $(cf["a_s_max"]):.3e}]===\n")
        my_lock = ReentrantLock()
        # Parallelize the loop
        @threads for iv in 1:length(v_0_list)
            for ia in 1:length(a_s_list)
                v_0 = v_0_list[iv]
                a_s = a_s_list[ia]
                v_0_norm = v_0 / e_perp
                gamma = -a_s * cf["n_atoms"] / l_perp
                sim.g = gamma2g(gamma, sim)
                sim.sigma2 = init_sigma2(sim.g)
                println(f">{ia:4d}, {iv:4d} | gamma = {gamma:.3e}, g = {sim.g:.3e}")
                sim.V0 = optical_lattice(v_0_norm,
                                          cf["d"] / l_perp,
                                          0.0,
                                          0.0,
                                          sim.X[1])
                if iv == 1 && ia == 1
                    sim.psi_0 = kspace(complex(gaussian(x / 3, sim)), sim)
                end
                sim.iswitch = -im

                # GS finding
                sol = runsim(sim, info=true)
                psi2 = zeros(cf["n"])
                try
                    sim.psi_0 = sol.u
                    psi2 = abs2.(xspace(sol.u, sim))
                    println("Success for $(ia), $(iv)")
                catch
                    println("Failure at $(ia), $(iv). Using initial guess.")
                    sim.psi_0 = kspace(complex(gaussian(x / 3, sim)), sim)
                end

                # Format output into a string and send to the channel
                output_str = "$ia, $iv | $(join(psi2, ", "))\n"
                @lock my_lock begin
                  write(file, f" {ia:4d}, {iv:4d} |")
                  for i in psi2
                    write(file, f" {i:10.9f}")
                  end
                  write(file, f"\n")
                  flush(file)
                end
            end
          end
    end
end

parallel_run()