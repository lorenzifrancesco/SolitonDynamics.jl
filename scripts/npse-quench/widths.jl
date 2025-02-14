#=
Script to reproduce the experiment for the width:
  - Obtian GS of pre-quench configuration
  - Quench and observe the result after 150ms
  - Compute the width of the result

Expected results:
  - region of collapse and explosion
  - region of stable shrinking

Input:   widths.csv, config.toml, config_pre_quench.toml
Output:  widths_experiment*.csv
=#
using FStrings
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


function estimate_width_rough(z_axis, final_psi2, dL)
  dz = z_axis[2] - z_axis[1]
  n_boxes = Int(round((z_axis[end]-z_axis[1]) / dL))
  box_z = LinRange(z_axis[1], z_axis[end], n_boxes)
  prob_box = Array{Float64, 1}(undef, n_boxes)
  for i in 1:n_boxes
    prob_box[i] = sum(final_psi2[(i-1)*dL .< z_axis .< i*dL]) * dz
    @info "box no. " i 
    @info "prob_box = " prob_box[i]
    @info "position = " box_z[i]

  end
  std = sqrt(sum(abs.(box_z) .* prob_box) / sum(prob_box))
  # print(f"\n center = {center:3.2e}, std = {std:3.2e} l_perp\n")
  return std
end


function estimate_width(z_axis, final_psi2)
  dz = z_axis[2] - z_axis[1]
  particle_fraction = sum(final_psi2) * dz
  # print(particle_fraction) # we can lose some particles
  center = sum(z_axis .* final_psi2) / particle_fraction
  std = sqrt(dz * sum(z_axis .^2 .* final_psi2) - sum(z_axis .* final_psi2)^2 / particle_fraction)
  print(f"\n center = {center:3.2e}, std = {std:3.2e} l_perp\n")
  return std
end


function particle_fraction(z_axis, final_psi2)
  dz = z_axis[2] - z_axis[1]
  return sum(final_psi2) * dz
end


# function get_widths()
begin
  #=
  Compute GS and post-quench dynamics for different values of a_s
  Save the final wavefunction for the state after 150ms. 
  Compute the width of the wavefunction.
  =#
  # constants: all units are in SI
  hbar_nostro = 1.0546e-34
  a_0 = 5.292e-11
  e_r = 2.1798723611030e-18

  info = false
  data_widths = CSV.read("input/widths.csv", DataFrame, header=false)
  rename!(data_widths, [:a_s, :width])
  # print(data_widths.a_s)

  cf_pre_quench = TOML.parsefile("input/config_pre_quench.toml")
  cf            = TOML.parsefile("input/config.toml")
  ### enforce structural properties of the normalization
  @assert(cf_pre_quench["omega_perp"]==cf["omega_perp"])
  @assert(cf_pre_quench["m"]==cf["m"])
  @assert(cf_pre_quench["n"]==cf["n"])
  @assert(cf_pre_quench["l"]==cf["l"])
  l_perp = sqrt(hbar_nostro/(cf_pre_quench["omega_perp"]*cf_pre_quench["m"]))
  e_perp = hbar_nostro * cf_pre_quench["omega_perp"]
  t_perp = cf_pre_quench["omega_perp"]^(-1) * 2 * pi # ATTENTION
  N = cf["n"]
  L = cf["l"]
  x = LinRange(-L / 2, L / 2, N + 1)[1:end-1]
  dx = x[2] - x[1]
  sim = init_sim((L,), (N,))

  params = data_widths.a_s
  result_widths = zeros(length(params))
  result_widths_rough = zeros(length(params))
  remaining_particle_fraction = zeros(length(params))

  ### Setup of the GS initialization pre-quench2
  # Note the parameters pre-quench2 remain fixed for
  # each post-quench configuration
  # g_maybe = 2 * hbar_nostro * cf_pre_quench["a_s"] * a_0 * cf_pre_quench["n_atoms"] / (cf_pre_quench["m"] * l_perp^2) # somehow not working TODO
  g_maybe = 2 * cf_pre_quench["a_s"] * a_0 * cf_pre_quench["n_atoms"] / (l_perp) # somehow not working TODO
  gamma = - cf_pre_quench["n_atoms"] * cf_pre_quench["a_s"] * a_0 / l_perp
  g = gamma2g(gamma, sim)
  @assert isapprox(g, g_maybe, atol=1e-9)
  sim.dt=1e-3
  sim.g = g
  sim.g5 = 0.0
  info && print("\n ------------------- GS pre-quench -------------------\n")
  l_x = sqrt(hbar_nostro/(cf_pre_quench["omega_x"]*cf_pre_quench["m"])) # SI
  e_recoil = (pi * hbar / cf_pre_quench["d"])^2 / (2*cf_pre_quench["m"]) # ATTENTION
  v_0_norm = cf_pre_quench["v_0"] * e_recoil / e_perp
  dL = cf_pre_quench["d"]
  print(f"\nl_perp = {l_perp *1e6 :6.3e} microns | dL = {dL* 1e6 :6.3e} microns = {dL/l_perp:.2e} L_perp\n")

  sim.V0 = optical_lattice(v_0_norm, 
                           cf_pre_quench["d"]/l_perp, 
                           0.0, 
                           l_x/l_perp,
                           sim.X[1])
  sim.sigma2 = init_sigma2(sim.g)
  sim.reltol = 1e-9
  sim.abstol = 1e-9
  sim.psi_0 = kspace(complex(gaussian(x/(3*cf_pre_quench["d"]/l_perp), sim)), sim)
  sim.iswitch = -im; 
  
  #### GS finding
  sim.maxiters = 1e6
  sol_gs = runsim(sim, info=false)
  
  # save for plotting
  psi2 = Array{Float64,2}(undef, (2, length(sol_gs.u)))
  for ix in [1, 2]
    psi2[ix, :] = abs2.(xspace(sol_gs.u, sim))
  end
  output = [vcat(sim.t[i], psi2[i, :]) for i in 1:size(psi2, 1)]
  CSV.write("results/widths_experiment_99.csv", DataFrame(output, :auto))

  print("\n")
  print("\n ------------------- Post-quench dynamics -------------------\n")
  for (idx, par) in enumerate(params)
    #### Computation of dynamics after the quench2
    sim.dt = 1e-3
    sim.psi_0 = sol_gs.u;
    gamma_post = - cf["n_atoms"] * par * a_0 / l_perp
    sim.g = gamma2g(gamma_post, sim)
    # @info "g = " g
    sim.sigma2 = init_sigma2(sim.g)
    e_recoil = (pi * hbar / cf["d"])^2 / cf["m"]
    v_0_norm = cf["v_0"] * e_recoil / e_perp
    l_x = sqrt(hbar_nostro/(cf["omega_x"]*cf["m"])) # SI
    # @info "v_0 norm =  " v_0_norm
    sim.V0 = optical_lattice(v_0_norm,
                             cf["d"]/l_perp, 
                             cf["tilt"], # tilt
                             l_x/l_perp,
                             sim.X[1])
    sim.iswitch = 1;
    sim.g5 = cf["l_3"] / l_perp^6 * t_perp * cf["n_atoms"]^3 / (6*pi^2)
    # @info "g5, 1D" sim.g5 
    sim.Nt = cf["n_t"]
    sim.tf = cf["t_f"] / t_perp
    # @info sim.tf
    sim.t = LinRange(sim.ti, sim.tf, sim.Nt)
    sim.time_steps = Int64(ceil((sim.tf-sim.ti)/sim.dt))
    sol = runsim(sim, info=false)

    info && print("\n ------------------- Saving -------------------\n")
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
    CSV.write("results/widths_experiment" * string(idx) * ".csv", DataFrame(output, :auto))
    
    if psi2[end, :] != zeros(N)
      result_widths[idx]       = estimate_width(real(sim.X[1]), psi2[end, :]) / (cf["d"]/l_perp)
      result_widths_rough[idx] = estimate_width_rough(real(sim.X[1]), psi2[end, :], cf["d"]/l_perp) / (cf["d"]/l_perp)
      remaining_particle_fraction[idx] = particle_fraction(real(sim.X[1]), psi2[end, :])
    else
      result_widths[idx]               = NaN
      result_widths_rough[idx]               = NaN
      remaining_particle_fraction[idx] = 0.0
    end
    print(f"({idx:5d} - {length(params):5d}) | as = {par:10.3e} a0 | width = {result_widths[idx]:10.3e} dL |  width_rough = {result_widths_rough[idx]:10.3e} dL | final_particles = {particle_fraction(real(sim.X[1]), psi2[end, :]):10.3e}\n")
  end
  CSV.write("results/widths_final.csv", DataFrame(a_s = params, width = data_widths.width, width_sim = result_widths, width_rough = result_widths_rough, particle_fraction = remaining_particle_fraction))
  return
end