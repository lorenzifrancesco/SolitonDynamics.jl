# Particular solvers to be used manually

# ============== Manual SplitStep methods, improved with exp
unpack_selection(sim, fields...) = map(x -> getfield(sim, x), fields)

function nlin_manual!(psi, sim::Sim{1,Array{ComplexF64}}, t; ss_buffer=nothing, info=false)
  g, X, V0, dV, equation, sigma2, dt, iswitch, N = unpack_selection(sim, :g, :X, :V0, :dV, :equation, :sigma2, :dt, :iswitch, :N)
  order = 2
  dt_order = dt/order
  x = X[1]
  N = N[1]
  xspace!(psi, sim)
  # 1D-GPE
  if equation == GPE_1D
    @. psi *= exp(dt_order * -im * iswitch * (V0 + V(x, t) + g * abs2(psi)))
  
  # NPSE
  elseif equation == NPSE
    nonlinear = g * abs2.(psi) ./ sigma2.(psi) + (1 ./ (2 * sigma2.(psi)) + 1 / 2 * sigma2.(psi))
    @. psi = exp(dt_order * -im * iswitch * (V0 + V(x, t) + nonlinear)) * psi
  
  # NPSE+
  elseif equation == NPSE_plus
    # load past solution
    if isnothing(ss_buffer)
      ss = ones(N)
    else
      ss = ss_buffer
    end
    sigma2_plus = zeros(length(x))
    M = N[1]
    dxx = 2*dV
    psisq = abs2.(psi)
    try
      # Nonlinear Finite Difference routine
      # ===================================
      function sigma_loop!(ret,sigma, params)
        # structure: [NPSE] + [simple derivatives of sigma] + [derivatives involving psi^2]
        @inbounds for j in 2:M-1
          ret[j] = (- sigma[j] .^ 4 + (1 + g*psisq[j])) - ((sigma[j+1]-sigma[j-1])/dxx)^2 +  sigma[j] * ((sigma[j-1]-2*sigma[j]+sigma[j+1])/(dV^2)) + sigma[j]*(sigma[j+1]-sigma[j-1])/dxx * (psisq[j+1]-psisq[j-1])/(dxx*psisq[j]) + (sigma[j+1]-sigma[j-1])/dxx * sigma[j] * (psisq[j+1]-psisq[j-1])/(dxx*psisq[j])
        end
        ret[1] = (- sigma[1] .^ 4 + (1 + g*psisq[1])) + ((sigma[2]-1.0)/dxx)^2 + ((1.0-2*sigma[1]+sigma[2])/(dV^2)) * sigma[1] +  (sigma[2]-1.0)/dxx * sigma[1] * (psisq[2]-0.0)/(dxx*psisq[1])
        ret[M] = (- sigma[M] .^ 4 + (1 + g*psisq[M])) - ((1.0-sigma[M-1])/dxx)^2 + ((sigma[M-1]-2*sigma[M]+1.0)/(dV^2)) * sigma[M] + (1.0-sigma[M-1])/dxx * sigma[M] * (0.0- psisq[M-1])/(dxx*psisq[M])
      end
      prob = NonlinearSolve.NonlinearProblem(sigma_loop!, ss, 0.0)
      sol = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), reltol=1e-6)
      sigma2_plus = (sol.u) .^ 2
      ss_buffer .= sol.u
    catch err
      if isa(err, DomainError)
        sigma2_plus = NaN
        throw(NpseCollapse(-666))
      else
        throw(err)
      end
    end
    temp = copy(sigma2_plus)
    temp = sqrt.(temp)
    temp_diff = copy(temp)
    # generate symmetric difference 
    temp_diff[1] = (temp[2]-temp[1])/dV 
    temp_diff[M] = (temp[M]-temp[M-1])/dV
    for i in 2:M-1
      temp_diff[i] = (temp[i+1]-temp[i-1])/dxx
    end
    nonlinear = g * abs2.(psi) ./ sigma2_plus + (1 / 2 * sigma2_plus) .+ (1 ./ (2 * sigma2_plus)) .* (1 .+ (temp_diff.^ 2))
    @. psi = exp(dt_order * -im * iswitch * (V0 + V(x, t) + nonlinear)) * psi

  # CQGPE
  elseif equation == CQGPE
    @. psi *= exp(dt_order * -im * iswitch * (V0 + V(x, t) + g * abs2(psi) - 6*log(4/3) * g^2 * abs2(abs2(psi))))
  end

  if maximum(abs2.(psi) * dV) > 0.8
    throw(Gpe3DCollapse(maximum(abs2.(psi) * dV)))
  end
  kspace!(psi, sim)
  return nothing
end


function propagate_manual!(psi, sim::Sim{1,Array{ComplexF64}}, t; ss_buffer=nothing, info=false)
  (ksquared, iswitch, mu, gamma_damp, dt) = unpack_selection(sim, :ksquared, :iswitch, :mu, :gamma_damp, :dt)
  psi_i = copy(psi)
  # splitting: N/2, N/2, L
  @. psi = exp(dt * iswitch * (1.0 - im * gamma_damp) * (-im * (1 / 2 * ksquared - mu))) * psi
  nlin_manual!(psi, sim, t; ss_buffer=ss_buffer, info=info)
  nlin_manual!(psi, sim, t; ss_buffer=ss_buffer, info=info)
  if iswitch == -im
    psi .= psi / sqrt(nsk(psi, sim))
    info && print(" - schempot: ", abs(chempotk_simple(psi, sim)))
    cp_diff = (chempotk_simple(psi, sim) - chempotk_simple(psi_i, sim)) / chempotk_simple(psi_i, sim) / dt
    return cp_diff
  else
    return nothing
  end
end

# ============== Manual CN GS

"""
Imaginary time evolution in xspace, 
using Crank Nicholson standard scheme
"""
function cn_ground_state!(psi, sim::Sim{1,Array{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
  @unpack dt, g, X, V0, iswitch, dV, Vol = sim
  x = X[1]
  psi_i = copy(psi)
  nonlin = (dt / 2) * g * abs2.(psi)
  tri_fwd += Diagonal(nonlin) # TODO check nonlinearity here
  tri_bkw += Diagonal(-nonlin)
  psi .= tri_fwd * psi
  psi .= transpose(\(psi, tri_bkw))
  psi .= psi / sqrt(ns(psi, sim))
  cp_diff = (chempot_simple(psi, sim) - chempot_simple(psi_i, sim)) / chempot_simple(psi_i, sim) / dt
  return cp_diff
end

# ============== Manual PC GS

"""
Imaginary time evolution in xspace, 
using predictor-corrector scheme (i FWD Euler + 1 fix point iterate)
"""
function pc_ground_state!(psi, sim::Sim{1,Array{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
  @unpack dt, g, X, V0, iswitch, dV, Vol, N = sim
  x = X[1]
  psi_i = copy(psi)
  nonlin = -(dt / 2) * g * abs2.(psi)
  tri_fwd += -Diagonal(ones(N[1])) + Diagonal(nonlin)
  for i in 1:3
    mapslices(x -> \(x, tri_fwd[i]), psi, dims=(i)) # FIXME am I a 3d method?
  end
  psi_star = tri_fwd * psi + psi
  psi .= 1 / 2 * (tri_fwd * psi) + psi

  nonlin_1 = -(dt / 2) * g * abs2.(psi_star)
  tri_fwd .= Diagonal(nonlin_1 - nonlin)
  psi .+= 1 / 2 * (tri_fwd * psi_star)
  tri_fwd = -Diagonal(ones(N[1])) + Diagonal(nonlin)
  psi .= 1 / 2 * (tri_fwd * psi_i + tri_fwd * psi) + psi
  info && @info display(sum(psi))

  psi .= psi / sqrt(ns(psi, sim))
  cp_diff = (chempot_simple(psi, sim) - chempot_simple(psi_i, sim)) / chempot_simple(psi_i, sim) / dt
  return cp_diff
end

# ============== Manual BE GS

"""
Imaginary time evolution in xspace, 
using BKW Euler
"""
function be_ground_state!(psi, sim::Sim{1,Array{ComplexF64}}, dt, tri_fwd, tri_bkw; info=false)
  @unpack dt, g, X, V0, iswitch, dV, Vol, N = sim
  x = X[1]
  psi_i = copy(psi)
  nonlin = -dt * (g * abs2.(psi) + V0)
  tri_bkw_complete = tri_bkw - Diagonal(nonlin)
  psi .= transpose(\(psi, tri_bkw_complete))
  psi .= psi / sqrt(ns(psi, sim))
  cp_diff = (chempot_simple(psi, sim) - chempot_simple(psi_i, sim)) / chempot_simple(psi_i, sim) / dt
  return cp_diff
end
