# Particular solvers to be used manually

# ============== Manual SplitStep methods, improved with exp

function nlin_manual!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack ksquared,g,X,V0,dV,Vol,mu,equation,sigma2,dt,iswitch = sim; x = X[1]
   xspace!(psi,sim)
   if equation == GPE_1D
      @. psi = exp(dt/2 * -im*iswitch* (V0 + V(x, t) + g*abs2(psi))) * psi
   elseif equation == NPSE
      nonlinear = g*abs2.(dpsi) ./sigma2.(dpsi) + (1 ./(2*sigma2.(dpsi)) + 1/2*sigma2.(dpsi))
      @. psi = exp(dt/2*-im*iswitch* (V0 + V(x, t) + nonlinear)) * psi
   end
   kspace!(psi,sim)
   return nothing
end


function propagate_manual!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma,dt = sim
   nlin_manual!(dpsi,psi,sim,t)
   @. psi = exp(dt/2 * (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)))*psi 
   return nothing
end


# ============== Manual SplitStep methods for ground state, improved with exp (wait for merge with above methods)

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
      @. psi += -dt/2*(1/2 * ksquared) * psi
      xspace!(psi, sim)
   return nothing
end

# ============== Manual CN GS

"""
Imaginary time evolution in xspace, 
using Crank Nicholson standard scheme
"""
function cn_ground_state!(psi,sim::Sim{1, Array{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
   @unpack dt,g,X,V0,iswitch,dV,Vol = sim; x = X[1]
   psi_i = copy(psi) 
   nonlin = (dt/2) * g*abs2.(psi)
   tri_fwd += Diagonal(nonlin)
   tri_bkw += Diagonal(-nonlin)
   #tri_fwd *= 1/2
   #tri_bkw *= 1/2
   psi .= tri_fwd*psi
   psi .= transpose(\(psi, tri_bkw))

   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end

# ============== Manual PC GS

"""
Imaginary time evolution in xspace, 
using predictor-corrector scheme (i FWD Euler + 1 fix point iterate)
"""
function pc_ground_state!(psi,sim::Sim{1, Array{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
   @unpack dt,g,X,V0,iswitch,dV,Vol,N = sim; x = X[1]
   psi_i = copy(psi) 
   nonlin = -(dt/2) * g*abs2.(psi)
   tri_fwd += - Diagonal(ones(N[1])) + Diagonal(nonlin)
   for i in 1:3
      mapslices(x -> \(x, tri_fwd[i]), psi, dims=(i))
   end
   psi_star = tri_fwd*psi + psi
   psi .= 1/2*(tri_fwd*psi) + psi

   nonlin_1 = -(dt/2) * g*abs2.(psi_star)
   tri_fwd .= Diagonal(nonlin_1-nonlin)
   psi .+= 1/2*(tri_fwd*psi_star) 
   tri_fwd = -  Diagonal(ones(N[1])) + Diagonal(nonlin)
   psi .= 1/2*(tri_fwd*psi_i + tri_fwd*psi) + psi
   @info display(sum(psi))

   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end

# ============== Manual BE GS

"""
Imaginary time evolution in xspace, 
using BKW Euler
"""
function be_ground_state!(psi,sim::Sim{1, Array{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
   @unpack dt,g,X,V0,iswitch,dV,Vol,N = sim; x = X[1]
   psi_i = copy(psi) 
   nonlin = -(dt/2) * g*abs2.(psi)
   tri_bkw += Diagonal(nonlin)
   psi .= transpose(\(psi, tri_bkw))

   norm_diff = ns(psi - psi_i, sim)/dt
   psi .= psi / sqrt(ns(psi, sim))
   return norm_diff
end
