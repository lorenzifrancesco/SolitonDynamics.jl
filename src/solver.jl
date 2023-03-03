function nlin!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack g,X,V0,iswitch,dV = sim; x = X[1]
   dpsi .= psi
   xspace!(dpsi,sim)

   if iswitch == -im
      @. dpsi *= -(V0 + V(x,t) + g*abs2(psi))
      mu2 = 1/(norm_squared(psi, sim)) .* (-dV .* sum(abs2.(psi) .* (V0 .+ V(x, t) .+ g .* abs2.(psi) )))
      @. dpsi += mu2*psi
   else
      @. dpsi *= -im * (V0 + V(x, t) + g*abs2(dpsi)) 
   end
   kspace!(dpsi,sim)
   return nothing
end

"""
Propagation routine
"""
function propagate!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, equation, dV, Vol = sim
   if equation == GPE_1D
      if iswitch == -im
         nlin!(dpsi,psi,sim,t)
         @. dpsi = (-1/2*ksquared)*psi
         mu1 = 1/(norm_squared(xspace(psi, sim), sim) .* Vol) .* sum(-1/2 .*ksquared .*abs2.(psi))
         @. dpsi += mu1*psi
      else
         nlin!(dpsi,psi,sim,t)
         @. dpsi += (-im*1/2*ksquared)*psi
      end
   elseif equation == NPSE
      throw("Unimplemented")
   end
   return nothing
end


function runsim(sim; info=false)
   @unpack psi_0, dV, ti, tf, solver = sim
   info && @info norm_squared(psi_0, sim)
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