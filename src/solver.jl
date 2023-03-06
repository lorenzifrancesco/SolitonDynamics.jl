"""
Time evolution in xspace
"""
function nlin!(dpsi,psi,sim::Sim{3, CuArray{ComplexF64}},t)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol = sim; x = X[1]
   dpsi .= psi
   ifft!(dpsi)
   @. dpsi *= -im * (V0 + V(x, t) + g*abs2(dpsi)) 
   fft!(dpsi)
   return nothing
end


"""
Time evolution in kspace
"""
function propagate!(dpsi, psi, sim::Sim{3, CuArray{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol = sim
   nlin!(dpsi,psi,sim,t)
   @. dpsi += -im*(1/2*ksquared)*psi
   return nothing
end

"""
Time evolution in xspace
"""
function nlin!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol = sim; x = X[1]
   dpsi .= psi
   xspace!(dpsi,sim)
   mu = 0.0
   if iswitch == -im
      mu = chempotk(psi, sim) # wainting for a more efficient implementation
   end
   @. dpsi *= -im *iswitch* (V0 + V(x, t) + g*abs2(dpsi)) + mu
   kspace!(dpsi,sim)
   return nothing
end


"""
Time evolution in kspace
"""
function propagate!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, equation, dV, Vol = sim
   if equation == GPE_1D
      nlin!(dpsi,psi,sim,t)
      @. dpsi += -im*iswitch*(1/2*ksquared)*psi
      #@info "norm" nsk(psi, sim)

   elseif equation == NPSE
      throw("Unimplemented")
   end
   return nothing
end

"""
Imaginary time evolution in xspace,
including explicit normalization
"""
function ground_state_nlin!(psi,sim::Sim{1, Array{ComplexF64}}, dt; info=false)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol = sim; x = X[1]
   
   psi_i = copy(psi) 
   ground_state_evolve!(psi, sim, dt)
   mu = 0.0
   if iswitch == -im
      mu = chempot(psi, sim) # wainting for a more efficient implementation
   end
   @. psi += - dt/2 * (V0 + g*abs2(psi)) * psi + dt/2*mu * psi

   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end

"""
Imaginary time evolution in kspace
"""
function ground_state_evolve!(psi, sim::Sim{1, Array{ComplexF64}}, dt; info=false)
      @unpack ksquared,g,X,V0,iswitch,dV,Vol = sim; x = X[1]
      kspace!(psi, sim)
      @. psi += -dt/2*(1/2 *ksquared) * psi
      xspace!(psi, sim)
   return nothing
end

"""
Main solution routine
"""
function runsim(sim; info=false)
   @unpack psi_0, dV, ti, tf, solver, iswitch = sim
   info && @info ns(psi_0, sim)

   # due to normalization, ground state solution 
   # is computed with forward Euler
   boring = true
   if iswitch == -im 
      if boring == false
         # xspace!(psi_0, sim)
         if solver == SplitStep 

            ssalg = DynamicSS(Tsit5(); abstol = 3e-3, 
            reltol = sim.reltol,
            tspan = Inf)

            problem = ODEProblem(propagate!, psi_0, (ti, tf), sim)
            ss_problem = SteadyStateProblem(propagate!, psi_0, sim)

            function normal!(u,t,integrator)
               @info "norm" nsk(u, sim)
            end
            cb = FunctionCallingCallback(normal!;
            func_everystep=true,
            func_start = true,
            tdir=1)

            sim.nfiles ?
            (sol = solve(ss_problem,
                        alg=ssalg,
                        callback=cb,
                        dense=false,
                        maxiters=1e8,
                        progress=true)) :
            (sol = solve(ss_problem,
                        alg=ssalg,
                        callback=cb,
                        dense=false,
                        maxiters=1e8,
                        progress=true))
         elseif solver == CrankNicholson
            throw("Unimplemented")
         end
      else
         xspace!(psi_0, sim)
         if solver == SplitStep 
            abstol_diff = 1e-8
            norm_diff = 1
            dt = 0.001
            #for i in  1:10000
            while norm_diff > abstol_diff
               norm_diff = ground_state_nlin!(psi_0,sim,dt)
               @info norm_diff
               # if norm_diff > abstol_diff * 1e5
               #    @warn "too fast"
               #    dt = dt * 0.9
               # end
            end 
         elseif solver == CrankNicholson
            throw("Unimplemented") 
         end
         kspace!(psi_0, sim)
         sol = psi_0
      end
      return [sol]
   else
      if solver == SplitStep 
         problem = ODEProblem(propagate!, psi_0, (ti, tf), sim)
         sim.nfiles ?
         (sol = solve(problem,
                     alg=sim.alg,
                     reltol=sim.reltol,
                     dense=false,
                     maxiters=1e10,
                     progress=true)) :
         (sol = solve(problem,
                     alg=sim.alg,
                     reltol=sim.reltol,
                     dense=false,
                     maxiters=1e10,
                     progress=true))

      elseif solver == CrankNicholson
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