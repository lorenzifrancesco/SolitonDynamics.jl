function nlin!(dphi,phi,sim::Sim{1, Array},t)
   @unpack g,X,V0 = sim; x = X[1]
   dphi .= phi
   xspace!(dphi,sim)
   @. dphi *= V0 + V(x,t) + g*abs2(dphi)
   kspace!(dphi,sim)
   return nothing
end

function nlin!(dphi,phi,sim::Sim{3, CuArray},t)
   @unpack g,X,V0 = sim; x,y,z = X
   y = y'; z = reshape(z,(1,1,length(z)))
   dphi .= phi
   xspace!(dphi,sim)
   @. dphi *= V0 + V(x,y,z,t) + g*abs2(dphi)
   kspace!(dphi,sim)
   return nothing
end

"""
Propagation routine
"""
function propagate!(dpsi, psi, sim::Sim{1}, t)
   @unpack gamma, equation = sim
   if equation == GPE_1D
      nlin!(dphi,phi,sim,t)
      @. dphi = -im*(1.0 - im*gamma)*dphi
   elseif equation == NPSE 
      throw("Unimplemented")
   end

   return nothing
end

function propagate!(dpsi, psi, sim::Sim{3}, t)
   @unpack gamma = sim
   nlin!(dphi,phi,sim,t)
   @. dphi = -im*(multiplier)*dphi
   return nothing
end

function runsim(sim; info=false)
   @unpack psi_0, tspan, reltol, solver = sim
   if solver == SplitStep 
      problem = ODEProblem(propagate!, psi_0, tspan)
      sol = solve(problem, Tsit5(), reltol=reltol)
   elseif solver == CrankNicholson
      throw("Unimplemented yet!!")
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