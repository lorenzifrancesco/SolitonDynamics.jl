
include("solvers_1D_manual.jl")
# include("solvers_3D_manual.jl")

function manual_run(
  sim;
  info=false,
  debug=false,
  throw_collapse=false,
  return_maximum=false)
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
    @assert solver == SplitStep
    cp_diff = 1
    abstol_diff = abstol
    cnt = 0
    info && print("\n")
    #
    decay = 0 * 1e-5
    debug && @info "setting exp decay rate to" decay
    # ss_buffer = ones(N[1])
    #
    minimum_evolution_time = 40.0
    #
    info && print("Iteration number")
    pr = Progress(minimum([maxiters, 5000]); dt=1)
    cp_diff = 1e300
    dump = Ref{Float64}(0.0)
    while cnt < maxiters &&
      (cnt * sim.dt < minimum_evolution_time || abs(cp_diff) > abstol_diff)
      try
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
      info && @printf("\riter = %1.3e - chempot_diff = %3.2e", cnt, cp_diff)
      catch err
          if isa(err, NpseCollapse) && !throw_collapse
              showerror(stdout, err)
          elseif isa(err, Gpe3DCollapse) && !throw_collapse
              showerror(stdout, err)
          else
              throw(err)
          end
          return nothing
      end
      cnt += 1
      update!(pr, cnt)
    end
    info && print("\n")
    # info && @info "Computation ended after iterations" cnt
    if cnt == sim.maxiters
      info && @warn "Maxiters reached"
    end
    # print(abs2.(xspace(psi, sim)))
    # print(sim.sigma2.(xspace(psi, sim)))
    sol = CustomSolution(u=psi, sigma=sim.sigma2.(xspace(psi, sim)),  t=t, cnt=cnt)
    return sol

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
    collapse_detected = false
    for i = 1:time_steps
      if !collapse_detected
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
            @info "Saving up to the collapse..."
            collapse_detected = true
          else
            throw(err)
          end
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
      else
        if t[save_counter] < solve_time_axis[i]
          if length(N) == 1
            collection[:, save_counter] = complex(zeros(length(psi)))
            collection_sigma[:, save_counter] = complex(ones(length(psi)))
          else
            throw("Not implemented")
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
    end

    # save the last item
    if !collapse_detected
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
    else
      if length(N) == 1
        collection[:, save_counter] = complex(zeros(length(psi)))
        collection_sigma[:, save_counter] = complex(ones(length(psi)))
      else
        throw("Not implemented")
        collection[:, :, :, save_counter] = psi
      end
      save_counter += 1
      time += dt
      if info
        cnt += 1
        update!(pr, cnt)
      end
    end

    # generate the solution
    if  length(N) == 1  
      sol =
          CustomSolution(
            u=[collection[:, k] for k = 1:Nt],
            sigma=[collection_sigma[:, k] for k = 1:Nt],
            t=t,
            cnt=time_steps)
      else
        sol =
          CustomSolution(
            u=[collection[:, :, :, k] for k = 1:Nt],
            t=t,
            cnt=time_steps)
    end
    return sol, max_prob * dV
  end
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
  @assert isapprox(nsk(psi_0, sim), 1.0, rtol=1e-9)
  @assert manual == true
  if return_maximum
    sol, max_prob = manual_run(sim; info=info, return_maximum=true)
    return sol, max_prob
  else
    sol = manual_run(sim; info=info, return_maximum=false)
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
