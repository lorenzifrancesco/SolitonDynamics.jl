"""
estimate the Gaussian sigma2 parameter from 3D data
"""

function estimate_sigma2(psi_k, sim::Sim{1,Array{ComplexF64}})
  @unpack equation, N, g, dV = sim
  psi = xspace(psi_k, sim)
  s2 = ones(N[1])
  if equation == GPE_1D
    return s2
  elseif equation == NPSE
    @. s2 = sqrt(1 + sim.g * abs2(psi))
  elseif equation == NPSE_plus
    # Nonlinear Finite Element routine
    b = (1 .+ g * abs2.(psi))
    b[1] += 1.0 * 1 / (4 * dV)
    b[end] += 1.0 * 1 / (4 * dV)

    a = ones(length(b))
    A0 = 1 / (2 * dV) * SymTridiagonal(2 * a, -a)

    ss = ones(N[1])
    prob = NonlinearProblem(sigma_eq, ss, [b, A0, dV])
    sol = solve(prob, NewtonRaphson(), reltol=1e-6)
    s2 = (sol.u) .^ 2
  end
  return s2
end

function sigma_eq(sigma, params)
  b = params[1]
  A0 = params[2]
  dV = params[3]
  N = length(sigma)
  bc = zeros(N)
  bc[1] = 1
  bc[end] = 1
  f = -1 / 2 * (A0 * sigma .^ 2) + 2 * sigma .* (A0 * sigma) + sigma .^ 4 - b - 1 / (2 * dV) * bc .* sigma
  return f
end

function estimate_sigma2(psi_k, sim::Sim{3,CuArray{ComplexF64}})
  s2 = Array{Float64,1}(undef, sim.N[1])
  # MSE estimator
  psi = xspace(psi_k, sim)
  dx = sim.X[1][2] - sim.X[1][1]
  aa = abs2.(psi)
  xax = 1:sim.N[1]
  yax = 1:sim.N[2]
  zax = 1:sim.N[3]
  tmp = zeros(sim.N[1])
  axial_density = sum(aa, dims=(2, 3))[:, 1, 1]
  ymask = (CuArray(sim.X[2]) .^ 2) * CuArray(ones(sim.N[2]))'
  zmask = CuArray(ones(sim.N[3])) * (CuArray(sim.X[3]) .^ 2)'
  r2mask = ymask+zmask
  display(sim.X[2])
  for x in xax
    if axial_density[x] < 1e-6
      tmp[x] = 1.0
      # @warn "found small prob"
    else
      tmp[x] = 1 - sum(aa[x, :, :] .* r2mask) * dx * dx # FIXME
    end
  end
  s2 = tmp
  return s2
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
compute normalization
"""
function ns(psi, sim)
  return sum(abs2.(psi)) * sim.dV
end

"""
compute normalization in k-space
"""
function nsk(psi, sim)
  dV, dk = measures(sim.L, sim.N)
  return sum(abs2.(psi)) * dk
end

"""
compute norm squared of a region
"""
function ns(psi, sim, mask)
  return sum(abs2.(psi) .* mask) * sim.dV
end

"""
return σ^2(ψ) of the NPSE
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
chemical potential of GPE in a given configuration
"""
function chempotk(psi, sim)
  @unpack ksquared, dV, V0, Vol, g = sim
  mu = 1 / Vol * sum(1 / 2 * ksquared .* abs2.(psi))
  tmp = xspace(psi, sim)
  mu += dV * sum((V0 + g * abs2.(tmp)) .* abs2.(tmp))
  mu *= 1 / ns(tmp, sim)
  #mu += 1 # add one transverse energy unit (1D-GPE case)
  return real(mu)
end

"""
chemical potential of GPE in a given configuration
"""
function chempot(psi, sim)
  @unpack ksquared, dV, V0, Vol, g = sim
  mu = dV * sum((V0 + g * abs2.(psi)) .* abs2.(psi))
  tmp = kspace(psi, sim)
  mu += 1 / Vol * sum(1 / 2 * ksquared .* abs2.(tmp))
  mu *= 1 / nsk(tmp, sim)
  #mu += 1 # add one transverse energy unit (1D-GPE case)
  return real(mu)
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
return the analytical (implicit) solution of the NPSE
"""
function npse_implicit(gamma, N)
  mu = npse_mu_full(gamma)
  positions = zeros(N)
  max = sqrt((1 - mu^2) / (2 * gamma))
  psi = LinRange(0, max, N) |> collect
  @. positions = 1 / sqrt(2) * 1 / sqrt(1 - mu) * atan(sqrt((sqrt(1 - 2 * gamma * abs2(psi)) - mu) / (1 - mu))) - 1 / sqrt(2) * 1 / sqrt(1 + mu) * atanh(sqrt((sqrt(1 - 2 * gamma * abs2(psi)) - mu) / (1 + mu)))
  return (positions, psi)
end

"""
in a given NPSE ground state, we expect to have a peak value of psi given by this
"""
function max_npse_psi2(gamma)
  mu = npse_mu_full(gamma)
  return (1 - mu^2) / (2 * gamma)
end