"""
estimate the Gaussian sigma2 parameter from 3D data
"""

function estimate_sigma2k(psi_k, sim::Sim{1,Array{ComplexF64}})
  @assert nsk(psi_k, sim) ≈ 1
  tmp = xspace(psi_k, sim)
  return estimate_sigma2(tmp, sim)
end

function estimate_sigma2(psi, sim::Sim{1,Array{ComplexF64}})
  @unpack equation, N, g, dV = sim
  @assert ns(psi, sim) ≈ 1
  s2 = ones(N[1])
  if equation == GPE_1D
    return s2
  elseif equation == NPSE
    try
      @. s2 = sqrt(1 + sim.g * abs2(psi))
    catch err
      if isa(err, DomainError)
        s2 = ones(length(psi)) * NaN
        throw(err)
        throw(NpseCollapse(sim.g * maximum(abs2.(psi))))
      end
    end
  elseif equation == NPSE_plus
    # Nonlinear Finite Element routine
    psisq = abs2.(psi)
    Nx = N[1]
    ss = ones(N[1])
    prob =
      NonlinearSolve.NonlinearProblem(sigma_loop_external!, ss, [psisq, dV, g])
    sol = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), reltol=1e-6)
    s2 = (sol.u) .^ 2
  end
  return s2
end

function estimate_sigma2k(psi_k, sim::Sim{3,CuArray{ComplexF64}})
  s2 = Array{Float64,1}(undef, sim.N[1])
  # MSE estimator
  psi = xspace(psi_k, sim)
  dx = real(sim.X[1][2] - sim.X[1][1])
  dy = real(sim.X[2][2] - sim.X[2][1])
  dz = real(sim.X[3][2] - sim.X[3][1])
  aa = abs2.(psi)
  xax = 1:sim.N[1]
  yaxis = real(sim.X[2])
  zaxis = real(sim.X[3])
  tmp::Array{Float64} = zeros(sim.N[1])
  axial_density = sum(aa, dims=(2, 3))[:, 1, 1]*dy*dz
  ymask = (CuArray(sim.X[2]) .^ 2) * CuArray(ones(sim.N[2]))'
  zmask = CuArray(ones(sim.N[3])) * (CuArray(sim.X[3]) .^ 2)'
  r2mask = ymask + zmask
  
  ## original
  # for x in xax
  #   if axial_density[x] < 1e-5
  #     tmp[x] = 1.0
  #     # @warn "found small prob"
  #   else
  #     tmp[x] = sum(aa[x, :, :] .* r2mask) / sum(aa[x, :, :])
  #   end
  # end
  # s2 = tmp

  # alternative
  for ix in xax  
    for (iy, y) in enumerate(yaxis)
      for (iz, z) in enumerate(zaxis)
          if axial_density[ix] < 1e-30
            tmp[ix] = 1.0
            # @warn "found small prob"
          end
        tmp[ix] += (y^2+z^2) * abs2(psi[ix, iy, iz]) / sum(aa[ix, :, :])
      end
    end
  end
  s2 = tmp *dy*dz
  return s2
end

function sigma_eq(ret, sigma, params)
  b = params[1]
  D1 = params[2]
  D2 = params[3]
  fterm = 0.0 * params[4]
  bc1 = params[5]
  bc2 = params[6]
  bc3 = params[7]
  d1sigma = D1 * sigma
  d2sigma = D2 * sigma
  # structure: NPSE + derivatives + boundary conditions
  ret .=
    (-sigma .^ 4 + b) +
    (-(d1sigma) .^ 2 + sigma .* d2sigma + sigma .* fterm .* d1sigma) +
    (bc2 - 2 * bc3 .* sigma + bc2 .* sigma + sigma .* fterm .* bc1)
end

function sigma_loop_external!(ret, sigma, params)
  psisq = params[1]
  dV = params[2]
  g = params[3]
  dxx = 2 * dV
  M = length(psisq)
  # structure: [NPSE] + [simple derivatives of sigma] + [derivatives involving psi^2]
  @inbounds for j = 2:M-1
    ret[j] =
      (-sigma[j] .^ 4 + (1 + g * psisq[j])) * psisq[j] -
      ((sigma[j+1] - sigma[j-1]) / dxx)^2 * psisq[j] +
      sigma[j] * ((sigma[j-1] - 2 * sigma[j] + sigma[j+1]) / (dV^2)) * psisq[j] +
      sigma[j] * (sigma[j+1] - sigma[j-1]) / dxx * (psisq[j+1] - psisq[j-1]) / (dxx)
  end
  ret[1] =
    (-sigma[1] .^ 4 + (1 + g * psisq[1])) * psisq[1]+
    ((sigma[2] - 1.0) / dxx)^2 * psisq[1] +
    ((1.0 - 2 * sigma[1] + sigma[2]) / (dV^2)) * sigma[1]*psisq[1] +
    (sigma[2] - 1.0) / dxx * sigma[1] * (psisq[2] - 0.0) / (dxx)
  ret[M] =
    (-sigma[M] .^ 4 + (1 + g * psisq[M])) * psisq[M]-
    ((1.0 - sigma[M-1]) / dxx)^2 * psisq[M] +
    ((sigma[M-1] - 2 * sigma[M] + 1.0) / (dV^2)) * sigma[M] * psisq[M] +
    (1.0 - sigma[M-1]) / dxx * sigma[M] * (0.0 - psisq[M-1]) / (dxx)
end

# sigma_eq_nf = NonlinearFunction(sigma_eq; jac=sigma_eq_jacobian)
# sigma_eq_fast = NonlinearFunction(sigma_eq; sparsity = )

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


function npse_expr(mu)
  f = ((1 - mu)^(3 / 2) - 3 / 2 * (1 - mu)^(1 / 2)) * (2 * sqrt(2)) / 3
  return f
end

"""
  Compute the ground state energy normalized to the harmonic energy unit
"""
function npse_energy(n, as)
  steps = 500
  dn = n / steps
  energy = 0
  for n_int in LinRange(0, n, steps)
    gamma = as * n_int
    a = roots(x -> npse_expr(x) + gamma, 0.5 .. 1, Newton, 1e-4)
    mu = mid(a[1].interval)
    energy += mu
  end
  return energy * dn
end

function gpe_energy(n, as)
  steps = 500
  dn = n / steps
  energy = 0
  for n_int in LinRange(0, n, steps)
    gamma = as * n_int
    mu = 1 - gamma^2 / 2
    energy += mu
  end
  return energy * dn
end

"""
  Compute the chemical potential
"""
function npse_mu(gamma)
  a = roots(x -> npse_expr(x) + gamma, 0.5 .. 1, Newton, 1e-4)
  mu = 0
  try
    mu = mid(a[1].interval)
  catch err
    throw(NpseCollapse)
  end
  return mu - 1
end

npse_mu_full(mu) = npse_mu(mu) + 1

function gpe_mu(n, as)
  return -(n * as)^2 / 8
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

"""
  Return σ^2(ψ) of the NPSE
"""
function init_sigma2(g::Float64)
  function sigma2(psi::ComplexF64)
    result = psi
    try
      result = sqrt(1 + g * abs2(psi))
    catch err
      if isa(err, DomainError)
        result = NaN
        throw(NpseCollapse(g * maximum(abs2.(psi))))
      else
        throw(err)
      end
    end
    return result
  end
  return sigma2
end

"""
  Chemical potential of GPE in a given configuration
"""
function chempotk(psi, sim)
  @assert nsk(psi, sim) ≈ 1
  tmp = xspace(psi, sim)
  return chempot(tmp, sim)
end

"""
  Chemical potential of GPE in a given configuration
"""
function chempot(psi, sim)
  @unpack ksquared, dV, V0, Vol, g, equation = sim
  @assert ns(psi, sim) ≈ 1
  if equation == GPE_1D || equation == GPE_3D
    mu = dV * sum((V0 + g * abs2.(psi)) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
    # mu *= 1 / nsk(tmp, sim)
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
    mat = Tridiagonal(-ones(length(s2) - 1), zeros(length(s2)), ones(length(s2) - 1))
    s2d = 1 / (2 * dV) * mat * s2
    s2d[1] -= 1 / (2 * dV) * s2[end]
    s2d[end] += 1 / (2 * dV) * s2[1]
    mu =
      dV * sum(
        (V0 + g * abs2.(psi) ./ s2 + 1 / 2 * (1 ./ s2 .* (1 .+ s2d .^ 2) + s2)) .*
        abs2.(psi),
      )
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

function sim_info(sim)
  @unpack_Sim N, L, manual = sim
  print("\nSimulation in $(length(N))D, with a number of steps $(N)\n")
end

function gpe_analytical(x, gamma; x0::Float64=0.0)
  @assert gamma > 0
  return sqrt(gamma / 2) * 2 / (exp(gamma * (x - x0)) + exp(-(x - x0) * gamma))
end

"""
  Return the analytical (implicit) solution of the NPSE
"""
function npse_implicit(gamma, N)
  mu = npse_mu_full(gamma)
  positions = zeros(N)
  max = sqrt((1 - mu^2) / (2 * gamma))
  psi = LinRange(0, max, N) |> collect
  @. positions =
    1 / sqrt(2) * 1 / sqrt(1 - mu) *
    atan(sqrt((sqrt(1 - 2 * gamma * abs2(psi)) - mu) / (1 - mu))) -
    1 / sqrt(2) * 1 / sqrt(1 + mu) *
    atanh(sqrt((sqrt(1 - 2 * gamma * abs2(psi)) - mu) / (1 + mu)))
  return (positions, psi)
end

"""
  In a given NPSE ground state, we expect to have a peak value of psi given by this
"""
function max_npse_psi2(gamma)
  mu = npse_mu_full(gamma)
  return (1 - mu^2) / (2 * gamma)
end
