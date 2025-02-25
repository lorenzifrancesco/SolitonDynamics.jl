using Printf
using TOML
using CSV, DataFrames, Tables, Printf, FFTW
using SolitonDynamics, Plots;
gr();


function optical_lattice(v_0, d, tilt, l_x, space)
  # See Eq.(3) of [PRA 75 033622 (2007)]
  # v_0 is in normalized units (hbar omega_perp)
  # notice that the k_L is different from the one of Stratclyde
  if l_x == 0.0
    return -v_0 * cos.(2 * pi / d * space) + tilt * space
  else
    return -v_0 * cos.(2 * pi / d * space) + tilt * space + 1/2 * (space/l_x).^2
  end
end

begin
  # constants: all units are in SI
  hbar_nostro = 1.0546e-34
  a_0 = 5.292e-11
  e_r = 2.1798723611030e-18
  cf_pre_quench = TOML.parsefile("input/config_pre_quench.toml")
  cf = TOML.parsefile("input/config.toml")
  
  ### enforce structural properties of the normalization
  @assert(cf_pre_quench["omega_perp"]==cf["omega_perp"])
  @assert(cf_pre_quench["m"]==cf["m"])
  @assert(cf_pre_quench["n"]==cf["n"])
  @assert(cf_pre_quench["l"]==cf["l"])
  l_perp = sqrt(hbar_nostro/(cf_pre_quench["omega_perp"]*cf_pre_quench["m"]))
  e_perp = hbar_nostro * cf_pre_quench["omega_perp"]
  t_perp = cf_pre_quench["omega_perp"]^(-1) * 2 * pi # ATTENTION 
  @info "l_perp" l_perp
  N = cf["n"]
  L = cf["l"]
  x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
  dx = x[2] - x[1]
  sim = init_sim((L,), (N,))

  params = [(0, 0)]

  ### Setup of the GS initialization pre-quench2
  # Note the parameters pre-quench2 remain fixed for
  # each post-quench configuration
  gamma = - cf_pre_quench["n_atoms"] * cf_pre_quench["a_s"] * a_0 / l_perp
  g = gamma2g(gamma, sim)
  
  y_values = []
  for (idx, par) in enumerate(params)
    sim.dt=1e-3
    sim.g = g
    sim.g5 = 0.0
    l_x = sqrt(hbar_nostro/(cf_pre_quench["omega_x"]*cf_pre_quench["m"])) # SI
    e_recoil = (pi * hbar / cf_pre_quench["d"])^2 / cf_pre_quench["m"]
    v_0_norm = cf_pre_quench["v_0"] * e_recoil / e_perp
    @info "v_0 norm =  " v_0_norm
    @info "gamma = " cf_pre_quench["a_s"] * a_0 * cf_pre_quench["n_atoms"] / l_perp
    sim.V0 = optical_lattice(v_0_norm, 
                             cf_pre_quench["d"]/l_perp, 
                             0.0, 
                             l_x/l_perp,
                             sim.X[1])
    sim.sigma2 = init_sigma2(sim.g)
    sim.reltol = 1e-9
    sim.abstol = 1e-9
    sim.psi_0 = kspace(complex(gaussian(x/5, sim)), sim)
    sim.iswitch = -im; 
    
    # GS finding
    sim.maxiters = 1e6
    @info "computing GS"
    sol = runsim(sim, info=true)
    
    ### Computation of dynamics after the quench2
    sim.dt = 1e-3
    sim.psi_0 = sol.u;
    gamma_post = - cf["n_atoms"] * cf["a_s"] * a_0 / l_perp
    sim.g = gamma2g(gamma_post, sim)
    @info "g = " g
    sim.sigma2 = init_sigma2(sim.g)
    e_recoil = (pi * hbar / cf_pre_quench["d"])^2 / (2*cf_pre_quench["m"]) # ATTENTION
    v_0_norm = cf["v_0"] * e_recoil / e_perp
    l_x = sqrt(hbar_nostro/(cf["omega_x"]*cf["m"])) # SI
    @info "v_0 norm =  " v_0_norm
    sim.V0 = optical_lattice(v_0_norm,
                             cf["d"]/l_perp, 
                             cf["tilt"], # tilt
                             l_x/l_perp,
                             sim.X[1])
    sim.iswitch = 1;
    sim.g5 = cf["l_3"] / l_perp^6 * t_perp * cf["n_atoms"]^3 / (6*pi^2)
    @info "g5, 1D" sim.g5 
    sim.Nt = cf["n_t"]
    sim.tf = cf["t_f"] / t_perp
    @info sim.tf
    sim.t = LinRange(sim.ti, sim.tf, sim.Nt)
    sim.time_steps = Int64(ceil((sim.tf-sim.ti)/sim.dt))
    sol = runsim(sim, info=true)
    # print(sol)

    print("\n ------------------- Plotting -------------------\n")
    psi2 = Array{Float64,2}(undef, (length(sol[1].u), length(sol[1].u[1])))
    sigma = Array{Float64,2}(undef, (length(sol[1].sigma), length(sol[1].sigma[1])))
    # px_psi = Array{Float64,2}(undef, (length(sol[1].u), length(sol[1].u[1])))
    for (ix, u_t) in enumerate(sol[1].u)
      # psi2[ix, :] = vcat(real(im * diff(xspace(u_t, sim)) / dx), 0.0)
      psi2[ix, :] = abs2.(xspace(u_t, sim))
      sigma[ix, :] = abs.(sim.sigma2.(xspace(u_t, sim)))
    end
    output = [vcat(sim.t[i], psi2[i, :]) for i in 1:size(psi2, 1)]
    output_sigma = [vcat(sim.t[i], sigma[i, :]) for i in 1:size(sigma, 1)]
    CSV.write("results/experiment" * string(idx) * ".csv", DataFrame(output, :auto))
    print(sum(sigma)/length(sigma))
  end
end