using ColorSchemes

# compute ground states
function check_1d_correctness(; gamma=0.0, info=true, iterate_t=false, show_waves=false, show=true)
  # TODO iterate here about N
  N_range = collect(200:200:20000)
  t_range = [0.05]
  mus = zeros(length(N_range), length(t_range))
  linf = zeros(length(N_range), length(t_range))
@warn mus
  for (iN, Nx) in enumerate(N_range)
    sd = load_parameters_alt(
      gamma_param=gamma,
      N_axial_1D=Nx,
      nosaves=true)
    sim = sd["G1"]

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
      analytical_sol = exp.(-x.^2 / 2) / sqrt(pi)
      p = plot(x, abs2.(analytical_sol), label="gaussian", color=:black)
    end
    analytical_sol = analytical_sol / sqrt(ns(analytical_sol, sim))

    # initialize psi_0 at random
    psi_0 .= complex(rand(N))
    psi_0 /= sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)

    @pack_Sim! sim

    plot!(p, x, abs2.(xspace(sim.psi_0, sim)), label="initial", color=:grey)
    nn = 2
    mm = 5
    pal = palette([:red, :blue], nn)
    for (idt, dt_set) in enumerate(t_range)
      @unpack_Sim sim
      dt = dt_set # XXX as in Triaxial
      abstol = 1e-6
      reltol = 1e-6
      @pack_Sim! sim

      sol = runsim(sim; info=info)
      plot!(p, x, abs2.(xspace(sol.u, sim)), color=pal[1], label="calculated")
      mus[iN, idt] = chempotk(sol.u, sim)        # TODO measure the final tau
      linf[iN, idt] = l_inf(sol.u - analytical_sol)
    end
    show_waves ? display(p) : nothing
  end

  if show
    p = plot(N_range, mus, label="mu", color=:red)
    display(p)
  end

  return (mus, linf)
end

function check_3d_correctness(; gamma=0.65, info=true, iterate_t=false, show_waves=false, show=true)
  # TODO iterate here about N
  N_3d = 512
  N_tran_range = collect(20:40:512)
  t_range = [0.05]
  mus = zeros(length(N_tran_range), length(t_range))
  linf = zeros(length(N_tran_range), length(t_range))

  for (iN, Nx) in enumerate(N_tran_range)
    sd = load_parameters_alt(
      gamma_param=gamma,
      N_axial_3D=N_3d,
      N_trans_3D=Nx,
      nosaves=true)
      
    sim = sd["G3"]

    @unpack_Sim sim
    x = X[1] |> real
    solver = SplitStep

    @assert length(sim.N) == 3
    @assert gamma > 0.0
    # ## get the analytical ground state
    # if gamma > 0.0
    #   analytical_sol = gpe_analytical.(x, gamma; x0=sim.L[1] / 4)
    #   true_min = chempot(gpe_analytical.(x, gamma; x0=sim.L[1] / 4), sim)
    #   p = plot(x, abs2.(gpe_analytical.(x, gamma; x0=sim.L[1] / 4)), label="soliton", color=:black)
    # else
    #   # harmonic oscillator
    #   @. V0 = 1 / 2 * (x^2)
    #   analytical_sol = exp.(-x.^2 / 2) / sqrt(pi)
    #   p = plot(x, abs2.(analytical_sol), label="gaussian", color=:black)
    # end
    # analytical_sol = analytical_sol / sqrt(ns(analytical_sol, sim))

    # initialize psi_0 at random
    psi_0 .= complex(rand(N))
    psi_0 /= sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)

    @pack_Sim! sim

    # plot!(p, x, abs2.(xspace(sim.psi_0, sim)), label="initial", color=:grey)
    nn = 2
    mm = 5
    pal = palette([:red, :blue], nn)
    for (idt, dt_set) in enumerate(t_range)
      @unpack_Sim sim
      dt = dt_set # XXX as in Triaxial
      abstol = 1e-6
      reltol = 1e-6
      @pack_Sim! sim

      sol = runsim(sim; info=info)
      # plot!(p, x, abs2.(xspace(sol.u, sim)), color=pal[1], label="calculated")
      mus[iN, idt] = chempotk(sol.u, sim)        # TODO measure the final tau
      linf[iN, idt] = l_inf(sol.u)
    end
    show_waves ? display(p) : nothing
    
    # lower the GC pressure
    sim = nothing
    sd = nothing
    sol = nothing
    
    GC.gc(true)
    CUDA.reclaim()
    CUDA.memory_status()
  end

  if show
    p = plot(N_tran_range, mus, label="mu", color=:red)
    display(p)
  end

  return (mus, linf)
end

function l_inf(x::AbstractArray)
  return maximum(abs.(x))
end