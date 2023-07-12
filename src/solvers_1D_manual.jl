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
  if equation == GPE_1D
    @. psi *= exp(dt_order * -im * iswitch * (V0 + V(x, t) + g * abs2(psi)))
  elseif equation == NPSE
    # @warn "whois sigma2? if appropriate this is zero" sigma2(Complex(sqrt(minimum(abs2.(psi)))))
    nonlinear = g * abs2.(psi) ./ sigma2.(psi) + (1 ./ (2 * sigma2.(psi)) + 1 / 2 * sigma2.(psi))
    @. psi = exp(dt_order * -im * iswitch * (V0 + V(x, t) + nonlinear)) * psi
  elseif equation == NPSE_plus
    sigma2_plus = zeros(length(x))
    try
      # Nonlinear Finite Element routine
      b = (1 .+ g * abs2.(psi))
      b[1] += 1.0 * 1 / (4 * dV)
      b[end] += 1.0 * 1 / (4 * dV)

      a = ones(length(b))
      A0 = 1 / (2 * dV) * SymTridiagonal(2 * a, -a)
      if isnothing(ss_buffer)
        ss = ones(N)
      else
        # info && @info "using ss_buffer with min: " minimum(ss_buffer)
        ss = ss_buffer
      end

      prob = NonlinearProblem(sigma_eq, ss, [b, A0, dV])
      sol = solve(prob, NewtonRaphson(), reltol=1e-6)
      sigma2_plus = (sol.u) .^ 2
      ss_buffer .= sol.u

      # === scientific debug zone
      append!(time_of_sigma, t)
      append!(sigma2_new, [sigma2_plus])
      append!(sigma2_old, [1 .+ g * abs2.(psi)])
      # === end scientific debug zone
    catch err
      if isa(err, DomainError)
        sigma2_plus = NaN
        throw(NpseCollapse(-666))
      else
        throw(err)
      end
    end
    temp = copy(sigma2_plus)
    nonlinear = g * abs2.(psi) ./ sigma2_plus + (1 / 2 * sigma2_plus .+ (1 ./ (2 * sigma2_plus)) .* (1 .+ (1 / dV * diff(prepend!(temp, 1.0))) .^ 2))
    @. psi = exp(dt_order * -im * iswitch * (V0 + V(x, t) + nonlinear)) * psi
    # Base.GC.gc()  
  elseif equation == CQGPE
    @. psi *= exp(dt_order * -im * iswitch * (V0 + V(x, t) + g * abs2(psi) - 1 / 4 * g^2 * abs2(abs2(psi))))
  end
  if maximum(abs2.(psi) * dV) > 0.5
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
    info && print(" - chempot: ", abs(chempotk(psi, sim)))
    cp_diff = (chempotk(psi, sim) - chempotk(psi_i, sim)) / chempotk(psi_i, sim) / dt
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
  cp_diff = (chempot(psi, sim) - chempot(psi_i, sim)) / chempot(psi_i, sim) / dt
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
    mapslices(x -> \(x, tri_fwd[i]), psi, dims=(i))
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
  cp_diff = (chempot(psi, sim) - chempot(psi_i, sim)) / chempot(psi_i, sim) / dt
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
  nonlin = -(dt / 2) * g * abs2.(psi)
  tri_bkw += Diagonal(nonlin)
  psi .= transpose(\(psi, tri_bkw))

  psi .= psi / sqrt(ns(psi, sim))
  cp_diff = (chempot(psi, sim) - chempot(psi_i, sim)) / chempot(psi_i, sim) / dt
  return cp_diff
end
