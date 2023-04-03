
# General solvers to be used with DiffEq


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
   elseif equation == NPSE_plus # compute the effect of the derivative perturbatively
      sigma2_0 = sigma2.(dpsi)
      # kspace!(sigma2_0, sim)
      # @. sigma2_0 *= -ksquared
      # xspace!(sigma2_0, sim)
      d2sigma = dV^2 * diff(diff(sigma2_0))
      append!(d2sigma, [0.0, 0.0])
      # test point: set the correction to zero and see NPSE collapse
      try
         sigma2_plus = sqrt.(-1/2 * d2sigma * 10  + g*abs2.(dpsi) .+ 1) # ad-hoc coefficient
      catch  err
         if isa(err, DomainError)
            sigma2_plus = NaN
            throw(NpseCollapse(g * maximum(abs2.(psi))))
         else
            throw(err)
         end
      end
      nonlinear = g*abs2.(dpsi) ./sigma2_plus + (1 ./(2*sigma2_plus) + 1/2*sigma2_plus)
      @. dpsi *= -im*iswitch* (V0 + V(x, t) + nonlinear) + mu_im
   end
   kspace!(dpsi,sim)
   return nothing
end


function propagate!(dpsi, psi, sim::Sim{1, Array{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma = sim
      nlin!(dpsi,psi,sim,t)
      #    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
      @. dpsi = (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)*psi + dpsi)
   return nothing
end