
"""
  Extract sigma^2 from a kspace wave function
"""
function estimate_sigma2k(psi_k, sim::Sim{1,Array{ComplexF64}})
  @assert nsk(psi_k, sim) ≈ 1
  tmp = xspace(psi_k, sim)
  return estimate_sigma2(tmp, sim)
end

"""
  Extract sigma^2 from a xspace wave function 
"""
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
  aa::Array{Float64} = Array(abs2.(psi))
  xax = 1:sim.N[1]
  yaxis = real(sim.X[2])
  zaxis = real(sim.X[3])
  tmp::Array{Float64} = zeros(sim.N[1])
  axial_density = sum(aa, dims=(2, 3))[:, 1, 1]*dy*dz
  @assert sum(axial_density*dx) ≈ 1


  ymask = (Array(sim.X[2]) .^ 2) * Array(ones(sim.N[2]))'
  zmask = Array(ones(sim.N[3])) * (Array(sim.X[3]) .^ 2)'
  r2mask = ymask + zmask
  

  #### original
  for x in xax
    if axial_density[x] < 1e-300
      tmp[x] = 1.0
      # @warn "found small prob"
    else
      tmp[x] = sum(aa[x, :, :] .* r2mask)*dy*dz / axial_density[x]
    end
  end
  s2 = tmp
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


"""
  Nonlinear problem corresponding to the NPSE+ second equation
"""
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
