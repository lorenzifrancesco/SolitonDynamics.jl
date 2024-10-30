
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

function gpe_mu(n, as)
  return -(n * as)^2 / 8
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
npse_mu_full(mu) = npse_mu(mu) + 1

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

function npse_expr(mu)
  f = ((1 - mu)^(3 / 2) - 3 / 2 * (1 - mu)^(1 / 2)) * (2 * sqrt(2)) / 3
  return f
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
