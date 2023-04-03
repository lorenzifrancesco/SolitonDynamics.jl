# Particular solvers to be used manually


"""
3D Imaginary time evolution in xspace, 
using BKW Euler
"""
function be_ground_state!(psi,sim::Sim{3, CuArray{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
   @unpack dt,g,X,V0,iswitch,dV,Vol,N = sim; x = X[1]
   throw("Broken")
   psi_i = copy(psi) 
   nonlin = -(dt/2) * g*abs2.(psi)
   tri_bkw += Diagonal(nonlin)
   psi .= transpose(\(psi, tri_bkw))

   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end