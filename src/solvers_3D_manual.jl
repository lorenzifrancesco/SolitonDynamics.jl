# Particular solvers to be used manually

# ============== Manual SplitStep methods, improved with exp

function nlin_manual!(psi,sim::Sim{3, CuArray{ComplexF64}},t; ss_buffer=nothing, info=false)
   @unpack ksquared,g,X,V0,dV,Vol,mu,equation,sigma2,dt,iswitch = sim; x = X[1]; y = X[1]; z = X[1]
   xspace!(psi,sim)
   @. psi *= exp(dt * -im*iswitch* (V0 + V(x,y,z,t)))
   @. psi *= exp(dt * -im*iswitch* (g*abs2(psi)))
   kspace!(psi,sim)
   return nothing
end

function propagate_manual!(psi, sim::Sim{3, CuArray{ComplexF64}}, t; ss_buffer=nothing, info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma_damp,dt = sim
   # info && @info "dt = " dt
   psi_i = copy(psi) 
   nlin_manual!(psi,sim,t; info=info)
   @. psi = exp(dt * iswitch * (1.0 - im*gamma_damp)*(-im*(1/2*ksquared - mu))) * psi
   if iswitch == -im
      psi .= psi / sqrt(nsk(psi, sim))
      info && print(" - chempot: ", chempotk(psi, sim))
      cp_diff = (chempotk(psi, sim) - chempotk(psi_i, sim))/(chempotk(psi_i, sim)) / dt
      return cp_diff
   else
      return nothing
   end
   return nothing
end

# ============== Manual BE GS

"""
3D Imaginary time evolution in xspace, 
using BKW Euler
"""
function be_ground_state!(psi,sim::Sim{3, CuArray{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
   @unpack dt,g,X,V0,iswitch,dV,Vol,N = sim; x = X[1]; y = X[1]; z = X[1]
   throw("Broken")
   psi_i = copy(psi) 
   nonlin = -(dt/2) * g*abs2.(psi)
   tri_bkw += Diagonal(nonlin)
   psi .= transpose(\(psi, tri_bkw))

   @warn "Still using old normalization"
   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end