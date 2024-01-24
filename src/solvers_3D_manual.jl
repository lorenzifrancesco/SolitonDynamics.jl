# Particular solvers to be used manually

# ============== Manual SplitStep methods, improved with exp

function nlin_manual!(
  psi::CuArray{ComplexF64},
  tmp_real1::CuArray{Float64},
  tmp_real2,
  sim::Sim{3,CuArray{ComplexF64}},
  t, 
  auxiliary;
  ss_buffer=nothing,
  info=false,
)

  @unpack ksquared, g, X, V0, dV, Vol, mu, equation, sigma2, dt, iswitch, collapse_threshold = sim
  order = 1
  dt_order = dt / order
  xspace!(psi, sim)
  @. psi *= exp(dt_order * -im * iswitch * (V0 + g * abs2(psi)))

  dydz = (X[2][2]-X[2][1])*(X[3][2]-X[3][1])
  tmp_real1 .= sum(abs2.(psi))*dydz
  if maximum(tmp_real1) > collapse_threshold / dV
    throw(Gpe3DCollapse(maximum(abs2.(psi) * dV)))
  end
  kspace!(psi, sim)
  return nothing
end

function propagate_manual!(
  psi::CuArray{ComplexF64},
  psi_i::CuArray{ComplexF64},
  tmp_real1::CuArray{Float64},
  tmp_real2::Array{Float64},
  sim::Sim{3,CuArray{ComplexF64}},
  t,
  aux;
  ss_buffer=nothing,
  info=false,
)
  @unpack ksquared, iswitch, dV, Vol, mu, gamma_damp, dt = sim
  @. psi =
    exp(dt * iswitch * (1.0 - im * gamma_damp) * (-im * (1 / 2 * ksquared - mu))) * psi
  nlin_manual!(psi, tmp_real1, tmp_real2, sim, t, aux; info=info)
  if iswitch == -im
    psi .= psi / sqrt(nsk(psi, sim))
    cp_diff =tmp_real1
      (chempotk_simple(psi, sim) - chempotk_simple(psi_i, sim)) /
      (chempotk_simple(psi_i, sim)) / dt
    psi_i .= psi
    return cp_diff
  else
    return 0.0
  end
  return 0.0
end

# ============== Manual BE GS

"""
3D Imaginary time evolution in xspace, 
using BKW Euler
"""
function be_ground_state!(
  psi,
  sim::Sim{3,CuArray{ComplexF64}},
  dt,
  tri_fwd,
  tri_bkw;
  info=false,
)
  @unpack dt, g, X, V0, iswitch, dV, Vol, N = sim
  x = X[1]
  y = X[1]
  z = X[1]
  throw("Broken")
  psi_i = copy(psi)
  nonlin = -(dt / 2) * g * abs2.(psi)
  tri_bkw += Diagonal(nonlin)
  psi .= transpose(\(psi, tri_bkw))

  @warn "Still using old normalization"
  norm_diff = ns(psi - psi_i, sim) / dt
  psi .= psi / sqrt(ns(psi, sim))
  return norm_diff
end
