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
  # TODO remove scalar indexing. Adapt the collapse 
  # dydz = (X[2][2] - X[2][1]) * (X[3][2] - X[3][1])
  tmp_real1 .= abs2.(psi) 
  # print(size(tmp_real1))
  # if maximum(tmp_real1) > collapse_threshold / dV
  #   throw(Gpe3DCollapse(maximum(abs2.(psi) * dV)))
  # end
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
    cp_diff = (chempotk_simple(psi, sim) - chempotk_simple(psi_i, sim)) /
              (chempotk_simple(psi_i, sim)) / dt
    psi_i .= psi
    return cp_diff
  else
    return 0.0
  end
end