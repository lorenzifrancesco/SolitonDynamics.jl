
include("solvers_1D_auto.jl")
include("solvers_1D_manual.jl")
include("solvers_3D_auto.jl")
include("solvers_3D_manual.jl")

function manual_run(
  sim;
  info=false,
  debug=false,
  throw_collapse=true,
  return_maximum=false,
)
  @unpack psi_0,
  dV,
  dt,
  ti,
  tf,
  t,
  solver,
  iswitch,
  abstol,
  reltol,
  N,
  Nt,
  V0,
  maxiters,
  time_steps,
  equation = sim
  psi = deepcopy(psi_0)

  tmp_psi = complex(zeros(N))
  tmp_real1 = zeros(N)
  if length(N)==3
    tmp_psi = CuArray(tmp_psi)
    tmp_real1 = CuArray(tmp_real1)
  end
  tmp_real2::Array{Float64} = zeros(N[1])
  sigma::Array{Float64} = ones(N[1])
  #######################
  # Imaginary time 
  #######################
  if iswitch == -im # select solver and run manual convergence routine 
    # in manual GS mode the maximum number of steps is specified by maxiters
    if solver == SplitStep
      cp_diff = 1
      abstol_diff = abstol
      cnt = 0
      info && print("\n")
      #
      debug && @info maxiters
      decay = 0 * 1e-5
      debug && @info "setting exp decay rate to" decay
      if equation == NPSE_plus
        ss_buffer = ones(N[1])
      else
        ss_buffer = nothing
      end
      #
      minimum_evolution_time = 40.0
      #
      info && print("Interaction number")
      pr = Progress(minimum([maxiters, 5000]); dt=1)
      cp_diff = 1e300
      dump = Ref{Float64}(0.0)
      while cnt < maxiters &&
        (cnt * sim.dt < minimum_evolution_time || abs(cp_diff) > abstol_diff)
        # try
        cp_diff = propagate_manual!(psi, 
                                    tmp_psi, 
                                    tmp_real1,
                                    tmp_real2, 
                                    sim, 
                                    0.0, 
                                    dump; 
                                    info=info, 
                                    ss_buffer=sigma)
        sim.dt *= (1 - decay)
        info && @printf("\riter = %5i - chempot_diff = %3.2e", cnt, cp_diff)
        # catch err
        #     if isa(err, NpseCollapse) && !throw_collapse
        #         showerror(stdout, err)
        #     elseif isa(err, Gpe3DCollapse) && !throw_collapse
        #         showerror(stdout, err)
        #     else
        #         throw(err)
        #     end
        #     return nothing
        # end
        cnt += 1
        update!(pr, cnt)
      end
      info && print("\n")
      info && @info "Computation ended after iterations" cnt
      sol = CustomSolution(u=[psi], t=t, cnt=cnt)
    else # nonspectral methods
      xspace!(psi, sim)
      solvers =
        [propagate_manual!, cn_ground_state!, pc_ground_state!, be_ground_state!]
      func = solvers[solver.number]
      info && @info "Solving using solver" func
      cp_diff = 1
      abstol_diff = abstol
      taglia = N[1]
      #for i in  1:10000
      d_central = -(dt / 2) * (1 / (dV^2) * ones(taglia)) |> complex
      d_lu = (dt / 4) * 1 / (dV^2) * ones(taglia - 1) |> complex
      tri_fwd = -SymTridiagonal(d_central, d_lu) + Diagonal(ones(taglia)) # Dx
      tri_bkw = SymTridiagonal(d_central, d_lu) + Diagonal(ones(taglia)) # Sx FIXME ??? sbagliato
      cnt = 0
      tmp = chempot_simple(psi, sim)
      while abs(cp_diff) > abstol_diff && cnt < maxiters
        cp_diff = func(psi, sim, dt, tri_fwd, tri_bkw)
        info && print("\n Interaction number")
        info && print("\r", cnt, " - chempot diff: ", cp_diff)
        # @assert tmp * cp_diff > 0
        tmp = cp_diff
        cnt += 1
      end
      kspace!(psi, sim)
      sol = CustomSolution(u=[psi], t=t, cnt=cnt)
    end
    return sol, -1.0

  #######################
  # Real time 
  ######################
  else
    # in manual run mode the number of steps is specified by time_steps
    time = 0.0
    max_prob = -1.0
    if length(N) == 1
      collection = Array{ComplexF64,2}(undef, (length(psi), Nt))
      collection = zeros((length(psi), Nt)) |> complex
      collection_sigma = Array{ComplexF64,2}(undef, (length(psi), Nt))
      collection[:, 1] = psi
      tmp_real1 = Array(zeros(N))
    else
      collection = CuArray{ComplexF64,4}(undef, (N..., Nt))
      collection[:, :, :, 1] = psi
      tmp_real1 = CuArray(zeros(N))
    end
    save_counter = 1
    solve_time_axis = LinRange(ti, tf, time_steps)
    #
    if equation == NPSE_plus
      ss_buffer = ones(N[1])
    else
      ss_buffer = nothing
    end
    debug && @warn "running with time_steps = " time_steps
    if info
      pr = Progress(time_steps)
      cnt = 0
    end
    auxiliary = 0.0
    auxiliary2 = Ref{Float64}(0.0)
    maximum_buffer = ones(N)
    for i = 1:time_steps
      try
        ## debug : feed ones after each repetition
        sigma .= ones(N[1])
        propagate_manual!(psi, 
                          tmp_psi, 
                          tmp_real1,
                          tmp_real2, 
                          sim, 
                          time, 
                          auxiliary2; 
                          info=info, 
                          ss_buffer=sigma)
        if return_maximum
          maximum_buffer = xspace(psi, sim)
          candidate_maximum = maximum(abs2.(maximum_buffer))
          if candidate_maximum > max_prob
            max_prob = candidate_maximum
          end
          if auxiliary < auxiliary2[]
            auxiliary = auxiliary2[]
            info && @info @sprintf("Maximum Linf error = %5.4e", auxiliary)
          end
        end
      catch err
        if isa(err, NpseCollapse) && !throw_collapse
          showerror(stdout, err)
        else
          throw(err)
        end
        return nothing
      end
      # print("\r", i, " - step")
      if t[save_counter] < solve_time_axis[i]
        if length(N) == 1
          collection[:, save_counter] = psi
          collection_sigma[:, save_counter] = sigma
        else
          collection[:, :, :, save_counter] = psi
        end
        save_counter += 1
      end
      time += dt
      if info
        cnt += 1
        update!(pr, cnt)
      end
    end
    if length(N) == 1
      collection[:, Nt] = psi
      collection_sigma[:, Nt] = sigma
      sol =
        CustomSolution(
          u=[collection[:, k] for k = 1:Nt],
          sigma=[collection_sigma[:, k] for k = 1:Nt],
          t=t,
          cnt=time_steps)
    else
      collection[:, :, :, Nt] = psi
      sol =
        CustomSolution(
          u=[collection[:, :, :, k] for k = 1:Nt],
          t=t,
          cnt=time_steps)
    end
    return sol, max_prob * dV
  end
end

function auto_run(sim; info=false, throw_collapse=true)
  @unpack psi_0,
  dV,
  dt,
  ti,
  tf,
  t,
  solver,
  iswitch,
  abstol,
  reltol,
  N,
  Nt,
  V0,
  maxiters,
  time_steps = sim
  psi = deepcopy(psi_0)
  @assert solver == SplitStep
  if iswitch == -im # solve a steady state problem
    #sim.iswitch = 1.0 # we should catch NPSE collapse in ground state?
    ssalg = DynamicSS(BS3(); reltol=sim.reltol, tspan=Inf)
    problem = ODEProblem(propagate!, psi, (ti, tf), sim)
    ss_problem = SteadyStateProblem(propagate!, psi, sim)
    sim.nfiles ?
    (
      sol = solve(
        ss_problem,
        alg=ssalg,
        callback=savecb,
        dense=false,
        maxiters=maxiters,
        progress=true,
        #dt = 0.001
      )
    ) :
    (
      sol = solve(
        ss_problem,
        alg=ssalg,
        dense=false,
        maxiters=maxiters,
        progress=true,
        #dt = 0.001
      )
    )
  else # propagate in real time
    problem = ODEProblem(propagate!, psi, (ti, tf), sim)
    try
      sim.nfiles ?
      (
        sol = solve(
          problem,
          alg=sim.alg,
          reltol=sim.reltol,
          saveat=sim.t[end],
          dt=dt,
          callback=savecb,
          dense=false,
          maxiters=maxiters,
          progress=true,
        )
      ) :
      (
        sol = solve(
          problem,
          alg=sim.alg,
          reltol=sim.reltol,
          saveat=sim.t,
          dt=dt,
          dense=false,
          maxiters=maxiters,
          progress=true,
        )
      )
    catch err
      if isa(err, NpseCollapse) && !throw_collapse
        showerror(stdout, err)
      else
        throw(err)
      end
      return nothing
    end
  end
  return sol
end

"""
Main solution routine
"""
function runsim(sim; info=false, return_maximum=false)
  @unpack psi_0,
  dV,
  dt,
  ti,
  tf,
  t,
  solver,
  iswitch,
  abstol,
  reltol,
  N,
  Nt,
  V0,
  maxiters,
  time_steps,
  manual = sim
  @assert !(!manual && return_maximum) # formal implication
  # function savefunction(psi...)
  #     isdir(path) || mkpath(path)
  #     i = findfirst(x -> x == psi[2], sim.t)
  #     padi = lpad(string(i), ndigits(length(sim.t)), "0")
  #     info && println("⭆ Save $i at t = $(trunc(ψ[2];digits=3))")
  #     # tofile = path*"/"*filename*padi*".jld2"
  #     tofile = joinpath(path, filename * padi * ".jld2")
  #     save(tofile, "ψ", psi[1], "t", psi[2])
  # end

  # savecb = FunctionCallingCallback(
  #     savefunction;
  #     funcat = sim.t, # times to save at
  #     func_everystep = false,
  #     func_start = true,
  #     tdir = 1,
  # )

  @assert isapprox(nsk(psi_0, sim), 1.0, rtol=1e-9)
  if manual == true
    if return_maximum
      sol, max_prob = manual_run(sim; info=info, return_maximum=true)
      return sol, max_prob
    else
      sol = manual_run(sim; info=info, return_maximum=false)
    end
  else
    sol = auto_run(sim; info=info)
  end
  return sol
end


function testsim(sim)
  err = false
  sol = try
    runsim(sim; info=true)
  catch e
    err = e
  end
  return sol, err
end
