using ColorSchemes

# compute ground states
function check_1d_space(;
  gamma=0.6,
  info=true,
  iterate_t=false,
  show_waves=false,
  show=true,
  eq="Np")

  N_range = [200, 400, 1024, 2048, 4096, 8192]
  dt_set = 0.05
  mus = zeros(length(N_range))
  elaps = zeros(length(N_range))
  linf = zeros(length(N_range))

  for (iN, Nx) in enumerate(N_range)
    sd = load_parameters_alt(
      gamma_param=gamma,
      N_axial_1D=Nx,
      nosaves=true)
    display(sd)
    sim = sd[eq]

    @unpack_Sim sim
    x = X[1] |> real
    solver = SplitStep
    manual = true
    @assert length(sim.N) == 1

    ## get the analytical ground state
    if gamma > 0.0
      analytical_sol = gpe_analytical.(x, gamma; x0=sim.L[1] / 4)
      true_min = chempot(gpe_analytical.(x, gamma; x0=sim.L[1] / 4), sim)
      p = plot(x, abs2.(gpe_analytical.(x, gamma; x0=sim.L[1] / 4)), label="soliton", color=:black)
    else
      # harmonic oscillator
      @. V0 = 1 / 2 * (x^2)
      analytical_sol = exp.(-x .^ 2 / 2) / sqrt(pi)
      p = plot(x, abs2.(analytical_sol), label="gaussian", color=:black)
    end
    analytical_sol = analytical_sol / sqrt(ns(analytical_sol, sim))

    # initialize psi_0 at random
    psi_0 .= complex(rand(N))
    psi_0 /= sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    dt = dt_set
    abstol = 1e-3
    reltol = 1e-3
    @pack_Sim! sim
    print("\n\n===> Relevant parameters: \n\t gamma = $gamma, N = $(sim.N), dt = $(sim.dt)")

    plot!(p, x, abs2.(xspace(sim.psi_0, sim)), label="initial", color=:grey)
    nn = 2
    pal = palette([:red, :blue], nn)

    elaps[iN] = @elapsed sol = runsim(sim; info=info)
    plot!(p, x, abs2.(xspace(sol.u, sim)), color=pal[1], label="calculated")
    mus[iN] = chempotk(sol.u, sim)
    linf[iN] = l_inf(sol.u - analytical_sol)
    show_waves ? display(p) : nothing
  end

  if show
    p = plot(N_range, mus, label="mu", color=:green, grid=:both)
    print("\n(dropping the first in elapsed time evaluation)")
    q = plot(N_range[2:end], elaps[2:end], label="execution_time", color=:red)
    display(q)
    display(p)
  end

  return (mus, linf)
end


function check_1d_time(;
  gamma=0.6,
  info=true,
  iterate_t=false,
  show_waves=false,
  show=true,
  eq="N")

  Nx = 2048
  dt_range = [0.2, 0.1, 0.05, 0.01, 0.001]
  mus = zeros(length(dt_range))
  elaps = zeros(length(dt_range))
  linf = zeros(length(dt_range))

  for (idt, dt_set) in enumerate(dt_range)
    sd = load_parameters_alt(
      gamma_param=gamma,
      N_axial_1D=Nx,
      nosaves=true)
    sim = sd[eq]

    @unpack_Sim sim
    x = X[1] |> real
    solver = SplitStep

    @assert length(sim.N) == 1

    ## get the analytical ground state
    if gamma > 0.0
      analytical_sol = gpe_analytical.(x, gamma; x0=sim.L[1] / 4)
      true_min = chempot(gpe_analytical.(x, gamma; x0=sim.L[1] / 4), sim)
      p = plot(x, abs2.(gpe_analytical.(x, gamma; x0=sim.L[1] / 4)), label="soliton", color=:black)
    else
      # harmonic oscillator
      @. V0 = 1 / 2 * (x^2)
      analytical_sol = exp.(-x .^ 2 / 2) / sqrt(pi)
      p = plot(x, abs2.(analytical_sol), label="gaussian", color=:black)
    end
    analytical_sol = analytical_sol / sqrt(ns(analytical_sol, sim))

    # initialize psi_0 at random
    psi_0 .= complex(rand(N))
    psi_0 /= sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    dt = dt_set
    abstol = 1e-3
    reltol = 1e-3
    @pack_Sim! sim

    print("\n\n===> Relevant parameters: \n\t gamma = $gamma, N = $(sim.N), dt = $(sim.dt)")
    plot!(p, x, abs2.(xspace(sim.psi_0, sim)), label="initial", color=:grey)
    nn = 2
    pal = palette([:red, :blue], nn)

    elaps[idt] = @elapsed sol = runsim(sim; info=info)
    plot!(p, x, abs2.(xspace(sol.u, sim)), color=pal[1], label="calculated")
    mus[idt] = chempotk(sol.u, sim)
    linf[idt] = l_inf(sol.u - analytical_sol)
    show_waves ? display(p) : nothing
    sim = nothing
  end

  if show
    p = plot(1 ./ dt_range, mus, label="mu", color=:green, xlabel="1/dt")
    display(p)
    print("\n(dropping the first in elapsed time evaluation)")
    q = plot(1 ./ dt_range, elaps, label="execution_time", color=:red, xlabel="1/dt")
    display(q)
  end

  return (mus, linf)
end

function check_3d_space(;
  gamma=0.65,
  info=true,
  show_waves=true,
  show=true)

  N_axial = 512
  N_tran_range = [60, 80, 120, 140]
  dt_set = 0.05
  mus = zeros(length(N_tran_range))
  linf = zeros(length(N_tran_range))

  p = plot()
  for (iN, Nx) in enumerate(N_tran_range)
    GC.gc(true)
    CUDA.reclaim()

    sd = load_parameters_alt(
      gamma_param=gamma,
      N_axial_3D=N_axial,
      N_trans_3D=Nx,
      nosaves=true,
      eqs=["G3"])
    sim = sd["G3"]
    print("\n============\n")
    CUDA.memory_status()
    print("\n============\n")

    @unpack_Sim sim
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    abstol = 1e-3
    reltol = 1e-3
    solver = SplitStep
    @assert length(sim.N) == 3
    @assert gamma > 0.0

    nn = 2
    mm = 5
    pal = palette([:red, :blue], nn)
    dt = dt_set
    @pack_Sim! sim

    print("\n\n===> Relevant parameters: \n\t gamma = $gamma, N = $(sim.N), dt = $(sim.dt)")
    
    sol = runsim(sim; info=info)

    @info "final imaginary time: " sim.dt * sol.cnt
    final = xspace(sol.u, sim)
    final_axial = Array(sum(abs2.(final), dims=(2, 3))[:, 1, 1] * (x[2] - x[1])^2)
    plot!(p, x, final_axial)

    mus[iN] = chempotk(sol.u, sim)
    linf[iN] = l_inf(sol.u)

    # lower the GC pressure
    sim = nothing
    sd = nothing
    sol = nothing

    GC.gc(true)
    CUDA.reclaim()
    CUDA.memory_status()
  end

  show_waves ? display(p) : nothing
  if show
    q = plot(N_tran_range, mus, label="mu", color=:red)
    display(q)
    savefig(q, "mu_vs_N_tran.pdf")
  end

  return (mus, linf)
end

function l_inf(x::AbstractArray)
  return maximum(abs.(x))
end