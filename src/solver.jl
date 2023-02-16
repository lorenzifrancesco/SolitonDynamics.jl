function nlin!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack g,X,V0,iswitch = sim; x = X[1]
   xspace!(dpsi,sim)
   psi .= dpsi
   if iswitch == -im
      @. dpsi *= iswitch*(V0 + V(x,t) + g*abs2(psi)) * psi + psi/(norm_squared(psi, sim)) * (-dV * sum(abs2.(psi) * (V(x, t) +g*abs2(psi) )))
   else
      @. dpsi *= (V0 + V(x,t) + g*abs2(psi)) * psi
   end
   kspace!(dpsi,sim)
   return nothing
end

"""
Propagation routine
"""
function propagate!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=true)
   @unpack ksquared, iswitch, equation, dV, V = sim
   if equation == GPE_1D
      if iswitch == -im
         nlin!(dpsi,psi,sim,t)
         @. dpsi = -1/2*ksquared*psi + psi/(norm_squared(xspace(psi, sim), sim) * V) * sum(-1/2*ksquared*psi) # normalization in reciprocal space warning!!
      else
         nlin!(dpsi,psi,sim,t)
         @. dpsi = -im*iswitch*ksquared*psi
      end
   elseif equation == NPSE
      throw("Unimplemented")
   end
   return nothing
end

function runsim(sim; info=true)
   @unpack psi_0, dV, ti, tf, solver = sim
   info && @info ob(psi_0, sim)

   if solver == SplitStep 
      problem = ODEProblem(propagate!, psi_0, (ti, tf), sim)
      sim.nfiles ?
      (sol = solve(problem,alg=sim.alg,saveat=sim.t[end],reltol=sim.reltol,callback=savecb,dense=false,maxiters=1e10,progress=true)) :
      (sol = solve(problem,alg=sim.alg,saveat=sim.t,reltol=sim.reltol,dense=false,maxiters=1e10,progress=true))
   elseif solver == CrankNicholson
      throw("Unimplemented")
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