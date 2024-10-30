"""
  Chemical potential in kspace
"""
function chempotk(psi, sim)
  @assert nsk(psi, sim) ≈ 1
  tmp = xspace(psi, sim)
  return chempot(tmp, sim)
end

"""
  Chemical potential in xspace
"""
function chempot(psi, sim)
  @unpack ksquared, dV, V0, Vol, g, equation = sim
  @assert ns(psi, sim) ≈ 1
  if equation == GPE_1D || equation == GPE_3D
    mu = dV * sum((V0 + g * abs2.(psi)) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
    if equation == GPE_1D
      mu += 1 # add one transverse energy unit (1D-GPE case)
    end
  elseif equation == NPSE
    s2 = estimate_sigma2(psi, sim)
    mu = dV * sum((V0 + g * abs2.(psi) ./ s2 + 1 / 2 * (1 ./ s2 + s2)) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
    # mu *= 1 / nsk(tmp, sim)
  elseif equation == NPSE_plus
    s2 = estimate_sigma2(psi, sim)
    M = length(psi)
    dxx = 2*dV
    tmp_real1 = copy(s2)
    s1 = sqrt.(s2) 
    tmp_real2 = copy(s1)
    tmp_real2[1] = (s1[2] - 1.0) / dxx
    tmp_real2[M] = (1.0 - s1[M-1]) / dxx
    @inbounds for i = 2:M-1
      tmp_real2[i] = (s1[i+1] - s1[i-1]) / dxx
    end
    nonlinear =
      g * abs2.(psi) ./ tmp_real1 + (1 / 2 * tmp_real1) .+
      (1 ./ (2 * tmp_real1)) .* (1 .+ (tmp_real2 .^ 2))
    mu = dV * sum((V0 + nonlinear) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
  elseif equation == CQGPE
    mu = dV * sum((V0 + g * abs2.(psi) - 1 / 4 * g^2 * abs2.(psi) .^ 2) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
    # mu *= 1 / nsk(tmp, sim)
    mu += 1 # add one transverse energy unit (1D-GPE case)
  end
  return real(mu)
end


"""
  Simplified (inaccurate) version of chempot,
  used in the ground state iterative step 
  Do not use for estimation of the physical chemical potential
"""
function chempotk_simple(psi, sim)
  @unpack ksquared, dV, V0, Vol, g, equation = sim
  mu::Float64 = 1 / Vol * sum(1 / 2 * ksquared .* abs2.(psi))
  tmp::AbstractArray{ComplexF64} = xspace(psi, sim)
  mu += dV * sum((V0 + g * abs2.(tmp)) .* abs2.(tmp))
  mu += 1 # add one transverse energy unit (1D-GPE case)
  return mu
end

function chempot_simple(psi, sim)
  @unpack ksquared, dV, V0, Vol, g, equation = sim
  mu::Float64 = dV * sum((V0 + g * abs2.(psi)) .* abs2.(psi))
  tmp::AbstractArray{ComplexF64} = kspace(psi, sim)
  mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
  mu += 1 # add one transverse energy unit (1D-GPE case)
  return mu
end

function project_radial(psi_k, sim::Sim{3,CuArray{ComplexF64}})
  # MSE estimator
  psi = xspace(psi_k, sim)
  aa = Array(abs2.(psi))
  x = sim.X[1] |> real
  y = sim.X[2] |> real
  z = sim.X[3] |> real
  rmax = sim.X[2][end] * sqrt(2) |> real
  r_steps = 128
  r_axis = LinRange(0, rmax, r_steps) |> collect
  dr = rmax / r_steps
  radial_density = zeros((sim.N[1], r_steps)) # axis and radius
  for (ix, xv) in enumerate(x)
    for (iy, yv) in enumerate(y)
      for (iz, zv) in enumerate(z)
        distance = sqrt(yv^2 + zv^2)
        ir = Int(round(r_steps * distance / rmax))
        radial_density[ix, ir] += aa[ix, iy, iz] * sim.dV
      end
    end
  end
  #radial_density[:, :] = aa[:, :, 64]
  radial_density .= radial_density / dr
  return r_axis, radial_density
end


"""
  Compute normalization
"""
function ns(psi, sim)
  return sum(abs2.(psi)) * sim.dV
end

"""
  Compute normalization in k-space
"""
function nsk(psi, sim)
  dV, dk = measures(sim.L, sim.N)
  return sum(abs2.(psi)) * dk
end

"""
  Compute norm squared of a region
"""
function ns(psi, sim, mask)
  return sum(abs2.(psi) .* mask) * sim.dV
end
