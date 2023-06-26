# General solvers to be used with DiffEq

function nlin!(dpsi,psi,sim::Sim{3, CuArray{ComplexF64}},t)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol,mu,equation,sigma2 = sim; x = X[1]; y = X[2]; z = X[3]
   dpsi .= psi
   mu_im = 0.0
   if iswitch == -im
      mu_im = chempotk(psi, sim) # wainting for a more efficient implementation
   end
   xspace!(dpsi,sim)
   @. dpsi *= -im*iswitch* (V0 + V(x,y,z,t) + g*abs2(dpsi)) + mu_im
   kspace!(dpsi,sim)
   return nothing
end


function propagate!(dpsi, psi, sim::Sim{3, CuArray{ComplexF64}}, t; info=false)
   @unpack ksquared, iswitch, dV, Vol,mu,gamma_damp = sim
      nlin!(dpsi,psi,sim,t)
      @. dpsi = (1.0 - im*gamma_damp)*(-im*(1/2*ksquared - mu)*psi + dpsi)
   return nothing
end