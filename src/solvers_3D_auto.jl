# General solvers to be used with DiffEq

#start old
# """
# Time evolution in xspace 3D
# """
# function nlin!(dpsi,psi,sim::Sim{3, CuArray{ComplexF64}},t)
#    @unpack g,X,V0,iswitch,dV,Vol,mu = sim; x = X[1]; y = X[2]; z = X[3]
#    dpsi .= psi
#    mu_im = 0.0
#    if iswitch == -im 
#       mu_im = chempotk(psi, sim) # wainting for a more efficient implementation
#    end
#    xspace!(dpsi, sim)
#    @. dpsi *= -im * (V0 + V.(x, y, z, t) + g*abs2(dpsi)) + mu_im
#    kspace!(dpsi, sim)
#    return nothing
# end


# """
# Time evolution in kspace 3D 
# """
# function propagate!(dpsi, psi, sim::Sim{3, CuArray{ComplexF64}}, t; info=false)
#    @unpack ksquared,iswitch,X,mu,gamma = sim; x = X[1]
#    nlin!(dpsi,psi,sim,t)
#    #@. dpsi = (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)*psi + dpsi)
#    @. dpsi += -im*(1/2*ksquared)*psi
#    #@info nsk(psi, sim)
#    return nothing
# end
#end

"""
Time evolution in xspace
"""
function nlin!(dpsi,psi,sim::Sim{3, CuArray{ComplexF64}},t)
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
      @. dpsi *= -im*iswitch* (V0 + V(x, y, z, t) + nonlinear) + mu_im
      #@info g*maximum(abs2.(dpsi))
   end
   kspace!(dpsi,sim)
   return nothing
end


"""
Time evolution in kspace
"""
function propagate!(dpsi, psi, sim::Sim{3, CuArray{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma = sim
      nlin!(dpsi,psi,sim,t)
      #    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
      @. dpsi = (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)*psi + dpsi)
   return nothing
end