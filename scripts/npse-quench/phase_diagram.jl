#=
Script to obtain all the NPSE ground states for the phase diagram.
  - Imaginary time solution for all the parameters.
  - Catch the NPSE / 3D GPE exceptions and save zeros.
  
Input:  phase_diagram.toml
Output: phase_diagram.txt (contains all the wavefunctions prob densities)
=#
using Printf
using TOML
using Base.Threads
using DelimitedFiles
using SolitonDynamics

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
  # constants: all units are in SI
  hbar = 1.0546e-34
  a_0 = 5.292e-11
  # e_r = 2.1798723611030e-18
  indices = []
  # GS phase diagram:
  cf = TOML.parsefile("input/phase_diagram.toml")
  e_recoil = (pi * hbar / cf["d"])^2 / (2 * cf["m"]) # ATTENTION using the notation of Strath
  l_perp = sqrt(hbar / (cf["omega_perp"] * cf["m"]))
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
  sim.iswitch = -im

  g_list   = LinRange(cf["g_min"] ,
                      cf["g_max"],
                      cf["n_g"])
  v_0_list = LinRange(cf["v_0_min"] * e_recoil / (hbar * cf["omega_perp"]), 
                      cf["v_0_max"] * e_recoil / (hbar * cf["omega_perp"]), 
                      cf["n_v_0"])
  v_0_min = cf["v_0_min"]
  v_0_max = cf["v_0_max"]
  g_min = cf["g_min"]
  g_max = cf["g_max"]
  print(v_0_list)
  open("results/phase_diagram.txt", "w") do file
    first_flag = true
    for (iv, v_0) in enumerate(v_0_list)
      for (ig, g) in enumerate(g_list)
        if !([ig, iv] in indices) 
          write(file, @sprintf(" %4d, %4d |", ig, iv))
          # gamma = - a_s * cf["n_atoms"] / (l_perp) #
          gamma = abs(g)/2 # ATTENTION, using the notaiton of Strath
          sim.g = gamma2g(gamma, sim)
          sim.sigma2 = init_sigma2(sim.g)
          @printf(">%4d, %4d | gamma = %.3e, g = %.3e\n", ig, iv, gamma, sim.g)
          sim.V0 = optical_lattice(v_0,
            cf["d"] / l_perp,
            0.0,
            0.0,
            sim.X[1])
          if first_flag == true
            print("RESETTING OF THE WF\n")
            sim.psi_0 = kspace(complex(gaussian(x/10, sim)), sim)
            first_flag = false
          end
          @warn( 1 + sim.g *  maximum(abs2.(xspace(sim.psi_0, sim))))
          @assert(1 + sim.g * maximum(abs2.(xspace(sim.psi_0, sim))) > 0)
          # print(sim)
          # GS finding
          psi2 = zeros(cf["n"])
          try
            sol = runsim(sim, info=false)
            sim.psi_0 = sol.u
            psi2 = abs2.(xspace(sim.psi_0, sim))
          catch e
            print("failed! ", e )
            first_flag = true
          end
          for i in psi2
            write(file, @sprintf(" %10.9f", i)) 
          end
          write(file, "\n")
          flush(file)
        else
          print("Skipping...")
        end
      end
    end
  end
end