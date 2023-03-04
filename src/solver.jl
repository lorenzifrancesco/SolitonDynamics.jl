"""
Time evolution in xspace
"""
function nlin!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol = sim; x = X[1]
   dpsi .= psi
   # compute the integral of mu for free
   mu_disp = sum(1/2 * ksquared .* abs2.(dpsi))
   xspace!(dpsi,sim)
   @. dpsi *= -im * (V0 + V(x, t) + g*abs2(dpsi)) 
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
      #    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)

      @. dpsi += -im*(1/2*ksquared)*psi
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
   
   initial_norm = ns(psi, sim) 
   ground_state_evolve!(psi, sim, dt)
   @. psi += - dt/2 * (V0 + g*abs2(psi)) * psi
   rel_change = abs(initial_norm - ns(psi, sim))/ns(psi, sim)
   psi .= psi / sqrt(ns(psi, sim))

   return rel_change
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
   if iswitch == -im
      xspace!(psi_0, sim)
      if solver == SplitStep 
         relative_tolerance = 6e-5
         rel_change = 1
         dt = 0.001
         #while rel_change > relative_tolerance
         for i in 1:10000
            rel_change = ground_state_nlin!(psi_0,sim,dt)
            if rel_change > 0.2
               @warn "too fast"
               dt = dt / 2
            end
         end 
      elseif solver == CrankNicholson
         throw("Unimplemented") 
      end
      kspace!(psi_0, sim)
      sol = [psi_0]
      return sol
   else
      if solver == SplitStep 
         problem = ODEProblem(propagate!, psi_0, (ti, tf), sim)
         sim.nfiles ?
         (sol = solve(problem,
                     alg=sim.alg,
                     reltol=sim.reltol,
                     callback=savecb,
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