# Particular solvers to be used manually

# ============== Manual SplitStep methods, improved with exp

function nlin_manual!(dpsi,psi,sim::Sim{3, CuArray{ComplexF64}},t)
   @unpack ksquared,g,X,V0,dV,Vol,mu,equation,sigma2,dt,iswitch = sim; x = X[1]; y = X[1]; z = X[1]
   xspace!(psi,sim)
   @. psi = exp(dt/2 * -im*iswitch* (V0 + V(x,y,z,t) + 2*pi*g*abs2(psi))) * psi
   kspace!(psi,sim)
   return nothing
end

function propagate_manual!(dpsi, psi, sim::Sim{3, CuArray{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma,dt = sim
   nlin_manual!(dpsi,psi,sim,t)
   @. psi = exp(dt/2 * (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)))*psi
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

   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end