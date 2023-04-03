
include("solvers_1D.jl")
include("solvers_1D ground_state.jl")
include("solvers_3D.jl")

"""
Main solution routine
"""
function runsim(sim; info=false)
   @unpack psi_0, dV, dt, ti, tf, t, solver, iswitch, abstol, reltol, N,Nt, V0, maxiters, time_steps, equation = sim
   info && @info ns(psi_0, sim)

   function savefunction(psi...)
      isdir(path) || mkpath(path)
      i = findfirst(x->x== psi[2],sim.t)
      padi = lpad(string(i),ndigits(length(sim.t)),"0")
      info && println("⭆ Save $i at t = $(trunc(ψ[2];digits=3))")
      # tofile = path*"/"*filename*padi*".jld2"
      tofile = joinpath(path,filename*padi*".jld2")
      save(tofile,"ψ",psi[1],"t",psi[2])
  end

  savecb = FunctionCallingCallback(savefunction;
                   funcat = sim.t, # times to save at
                   func_everystep=false,
                   func_start = true,
                   tdir=1)

   # due to normalization, ground state solution 
   # is computed with forward Euler
   boring = false
   if iswitch == -im 
      if boring == false
         sim.iswitch = 1.0
         if solver == SplitStep 

            ssalg = DynamicSS(BS3(); 
            reltol = sim.reltol,
            tspan = Inf)

            problem = ODEProblem(propagate!, psi_0, (ti, tf), sim)
            ss_problem = SteadyStateProblem(propagate!, psi_0, sim)

            sim.nfiles ?
            (sol = solve(ss_problem,
                        alg=ssalg,
                        callback=savecb,
                        dense=false,
                        maxiters=maxiters,
                        progress=true, 
                        #dt = 0.001
                        )) :
            (sol = solve(ss_problem,
                        alg=ssalg,
                        dense=false,
                        maxiters=maxiters,
                        progress=true, 
                        #dt = 0.001
                        ))
         elseif solver == CrankNicholson
            throw("Unimplemented")
         end
         return sol
      else
         xspace!(psi_0, sim)
         if solver == SplitStep 
            norm_diff = 1
            abstol_diff = abstol * dt
            #for i in  1:10000
            while norm_diff > abstol_diff
               norm_diff = ground_state_nlin!(psi_0,sim,dt)
               @info norm_diff
               # if norm_diff > abstol_diff * 1e5
               #    @warn "too fast"
               #    dt = dt * 0.9
               # end
            end
         else
            solvers = [ground_state_nlin!, cn_ground_state!, pc_ground_state!, be_ground_state!]
            func = solvers[solver.number]
            @info "Solving using solver" func 
            norm_diff = 1
            abstol_diff = abstol
            taglia = N[1]
            #for i in  1:10000
            d_central = -(dt/2) * ( 1/(dV^2) * ones(taglia) - V0) |> complex
            d_lu = (dt/2) * 1/(2*dV^2) * ones(taglia-1) |> complex
            tri_fwd = SymTridiagonal(d_central, d_lu) + Diagonal(ones(taglia)) # Dx
            tri_bkw = -SymTridiagonal(d_central, d_lu) + Diagonal(ones(taglia)) # Sx
            cnt = 0 
            while norm_diff > abstol_diff && cnt < maxiters
               norm_diff = func(psi_0,sim,dt, tri_fwd, tri_bkw)
               cnt +=1
            end
            @info "Computation ended after iterations" cnt
         end
         kspace!(psi_0, sim)
         sol = psi_0
         return [sol]
      end
   else # real-time dynamics
      if solver == SplitStep 
         if equation == NPSE_plus
            initial_sigma2 = 1.0
            initial_lambda = 0.0
            problem = ODEProblem(propagate!, psi_0, (ti, tf), sim, initial_sigma2, initial_lambda)
            try
            sim.nfiles ?
            (sol = solve(problem,
                        alg=sim.alg,
                        reltol=sim.reltol,
                        saveat=sim.t[end],
                        dt=dt,
                        callback=savecb,
                        dense=false,
                        maxiters=maxiters,
                        progress=true)) :
            (sol = solve(problem,
                        alg=sim.alg,
                        reltol=sim.reltol,
                        saveat=sim.t,
                        dt=dt,
                        dense=false,
                        maxiters=maxiters,
                        progress=true))
            catch err
               if isa(err, NpseCollapse)
                  showerror(stdout, err)
               else
                  throw(err)
               end
               return nothing
            end
         else
            problem = ODEProblem(propagate!, psi_0, (ti, tf), sim)
            try
            sim.nfiles ?
            (sol = solve(problem,
                        alg=sim.alg,
                        reltol=sim.reltol,
                        saveat=sim.t[end],
                        dt=dt,
                        callback=savecb,
                        dense=false,
                        maxiters=maxiters,
                        progress=true)) :
            (sol = solve(problem,
                        alg=sim.alg,
                        reltol=sim.reltol,
                        saveat=sim.t,
                        dt=dt,
                        dense=false,
                        maxiters=maxiters,
                        progress=true))
            catch err
               if isa(err, NpseCollapse)
                  showerror(stdout, err)
               else
                  throw(err)
               end
               return nothing
            end
         end
      elseif solver == ManualSplitStep
         time = 0.0
         psi = 0.0 * psi_0
         psi .= psi_0
         dpsi = 0.0 * psi
         collection = Array{ComplexF64, 2}(undef, (length(psi_0), Nt))
         collection[:, 1] = psi

         save_interval = Int(round(time_steps/Nt))
         for i in 1:time_steps
            propagate_manual!(dpsi, psi, sim, time)
            if i % save_interval == 0
               collection[:, Int(floor(i / save_interval))] = psi
            end
            @info "norm" (nsk(psi_0, sim))
            time += dt
         end
         sol = CustomSolution(u=[collection[:, k] for k in 1:Nt], t=t)

      elseif solver != SplitStep
         throw("Unimplemented")
      end
   end
   return sol
end


function testsim(sim)
   err = false
   sol = try
           runsim(sim; info=false)
       catch e
           err = true
       end
return sol,err
end