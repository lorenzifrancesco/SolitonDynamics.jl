
# General solvers to be used with DiffEq
"""
Time evolution in xspace
"""
function nlin!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol,mu,equation,sigma2 = sim; x = X[1]
   dpsi .= psi
   mu_im = 0.0
   if iswitch == -im
      mu_im = chempotk(psi, sim) # wainting for a more efficient implementation
   end
   xspace!(dpsi,sim)
   if equation == GPE_1D
      @. dpsi *= -im*iswitch* (V0 + V(x, t) + g*abs2(dpsi)) + mu_im
   elseif equation == NPSE
      nonlinear = g*abs2.(dpsi) ./sigma2.(dpsi) + (1 ./(2*sigma2.(dpsi)) + 1/2*sigma2.(dpsi))
      @. dpsi *= -im*iswitch* (V0 + V(x, t) + nonlinear) + mu_im
      #@info g*maximum(abs2.(dpsi))
   end
   kspace!(dpsi,sim)
   return nothing
end


"""
Time evolution in kspace
"""
function propagate!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma = sim
      nlin!(dpsi,psi,sim,t)
      #    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
      @. dpsi = (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)*psi + dpsi)
   return nothing
end


# ============== ManualSplitStep methods, improved with exp
"""
Time evolution in xspace
"""
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


"""
Time evolution in kspace
"""
function propagate_manual!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma,dt = sim
   nlin_manual!(dpsi,psi,sim,t)
   @. psi = exp(dt/2 * (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)))*psi 
   return nothing
end