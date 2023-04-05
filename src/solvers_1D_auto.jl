
# General solvers to be used with DiffEq

function nlin!(dpsi,psi,sim::Sim{1, Array{ComplexF64}},t)
   @unpack ksquared,g,X,V0,iswitch,dV,Vol,mu,equation,sigma2 = sim; x = X[1] |> real
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
      # test point: set the correction to zero and see NPSE collapse
      sigma2_plus = zeros(length(x))
      try
         init_sigma = 1.0
         init_lambda = 0.0
         problem = ODEProblem(sigma2_diff!, [init_sigma, init_lambda], (x[1], x[end]), [sim, dpsi])
         # initial state is (1, 0)
         sol = solve(problem,
                     alg=Tsit5(),
                     maxiters=1000,
                     saveat=x
                     )
         if length(sol.u) < length(dpsi)
            throw(DomainError("placeholder"))
            @warn "Danger length(sol.u) < length(dpsi)"
         end
         sigma2_plus = [sol.u[i][1] for i in 1:length(x)]

         # what is the difference in using this improved model?
         # try
         #    @info sum.(abs2.(sigma2_plus - sigma2.(psi))) # for each time step
         # catch err
         #    @info "collapse, pazienza"
         # end

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
      @. dpsi = (1.0 - im*gamma)*(-im*(1/2*ksquared - mu)*psi + dpsi)
   return nothing
end

function sigma2_diff!(dSigmaLambda, state, params, x)
   L = params[1].L[1]
   idx = minimum([Int(floor((x+L/2)/params[1].dV))+1, 256]) # horrible solution
   dSigmaLambda[1] = state[2]
   dSigmaLambda[2] = -2*state[1]^2 + 2*(1 + params[1].g * abs2(params[2][idx]))
   return nothing
end