# Particular solvers to be used manually

# ============== Manual SplitStep methods, improved with exp
unpack_selection(sim, fields...) = map(x -> getfield(sim, x), fields)

@inline function nlin_manual!(
  psi,
  tmp_real1,
  tmp_real2,
  sim::Sim{1,Array{ComplexF64}},
  t,
  auxiliary;
  ss_buffer=nothing,
  info=false,
)
  g, V0, dV, equation, sigma2, dt, iswitch, N, collapse_threshold, X =
    unpack_selection(sim, :g, :V0, :dV, :equation, :sigma2, :dt, :iswitch, :N, :collapse_threshold, :X)
  order = 1
  dt_order = dt / order
  N = N[1]
  xspace!(psi, sim)
  if equation == GPE_1D
    @. psi *= exp(dt_order * -im * iswitch * (V0 + g * abs2(psi)))
  elseif equation == NPSE
    nonlinear =
      g * abs2.(psi) ./ sigma2.(psi) +
      (1 ./ (2 * sigma2.(psi)) + 1 / 2 * sigma2.(psi))
    @. psi = exp(dt_order * -im * iswitch * (V0 + nonlinear)) * psi
  elseif equation == NPSE_plus
    M = N[1]
    dxx = 2 * dV
    # tmp_real1 <--- psi^2
    # tmp_real2 <--- derivative of psi^2
    tmp_real1 .= abs2.(psi)
    left_border = 1
    right_border = M
    cnt = 1
    if iswitch == 1 && false
      while tmp_real1[cnt] < 4e-5
        cnt += 1
      end
      left_border = cnt
      cnt = M
      while tmp_real1[cnt] < 4e-5
        cnt -= 1
      end
      right_border = cnt
    end
    # info && @info @sprintf("borders: %4i, %4i", left_border, right_border)

    ## Nonlinear Finite Difference routine
    ## ===================================
    try
      # ######################### BVProblem METHOD
      # # interpolation
      # tmp_real2[1] = (tmp_real1[2]) / dxx
      # tmp_real2[M] = (- tmp_real1[M-1]) / dxx
      # @inbounds for i = 2:M-1
      #   tmp_real2[i] = (tmp_real1[i+1] - tmp_real1[i-1]) / dxx
      # end
      # Xr = real.(X[1])
      # # info && @info Xr[1], Xr[end]
      # psisq_interp = Interpolations.linear_interpolation(Xr, tmp_real1)
      # psisq_derivative = Interpolations.linear_interpolation(Xr, tmp_real2)
      # function sigma_bvp!(du, u, p, x)
      #   # @info x
      #   du[1] = u[2]
      #   du[2] = u[1]^3 -
      #           (1+g*psisq_interp(x))/u[1] +
      #           u[2]^2/u[1] -
      #           u[2] * psisq_derivative(x)/psisq_interp(x)
      # end
      # function bc1!(residual, u, p, t)
      #     residual[1] = u[1][1] - 1.0
      #     residual[2] = u[end][1] - 1.0
      # end
      # bvp1 = BVProblem(sigma_bvp!, bc1!, [1.0, 0.0], [Xr[left_border], Xr[right_border]])
      # sol = solve(bvp1, BoundaryValueDiffEq.MIRK2(), dt=dV, adaptive=false, reltol=1e-3, abstol=1e-3)
      # # @info sol.retcode
      # ### check the correctness
      # for i in left_border:right_border
      #   ss_buffer[i] = sol.u[i+1-left_border][1]
      # end
      # # ################### END METHOD

      #################### Newton-Raphson METHOD
      ## check the variables
      psisq = tmp_real1[left_border:right_border]
      ## define function inside the restriction
      M_restr = length(psisq)
      function sigma_loop!(ret, sigma, params)
        # structure: [NPSE] + [simple derivatives of sigma] + [derivatives involving psi^2]
        @inbounds @simd for j = 2:M_restr-1
          ret[j] =
            (-sigma[j] .^ 4 + (1 + g * psisq[j])) + ((sigma[j+1]-sigma[j-1])/dxx)^2
            #* psisq[j] -
            # ((sigma[j+1] - sigma[j-1]) / dxx)^2 * psisq[j]+
            # sigma[j] * ((sigma[j-1] - 2 * sigma[j] + sigma[j+1]) / (dV^2)) * psisq[j]+
            # sigma[j] * (sigma[j+1] - sigma[j-1]) / dxx * (psisq[j+1] - psisq[j-1]) / (dxx)
        end
        ret[1] =
          (-sigma[1] .^ 4 + (1 + g * psisq[1])) + ((sigma[2]-1.0)/dxx)^2
          #* psisq[1] +
          # ((sigma[2] - 1.0) / dxx)^2*psisq[1] +
          # ((1.0 - 2 * sigma[1] + sigma[2]) / (dV^2)) * sigma[1]*psisq[1] +
          # (sigma[2] - 1.0) / dxx * sigma[1] * (psisq[2] - 0.0) / (dxx)
        ret[M_restr] =
          (-sigma[M_restr] .^ 4 + (1 + g * psisq[M_restr])) + ((1-0-sigma[M_restr-1])/dxx)^2
          #* psisq[M_restr] -
          # ((1.0 - sigma[M_restr-1]) / dxx)^2 *psisq[M_restr]+
          # ((sigma[M_restr-1] - 2 * sigma[M_restr] + 1.0) / (dV^2)) * sigma[M_restr] *psisq[M_restr]+
          # (1.0 - sigma[M_restr-1]) / dxx * sigma[M_restr] * (0.0 - psisq[M_restr-1]) / (dxx)
      end
      # jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> sigma_loop!(du, u, 0.0), ones(N[1]), ones(N[1]))
      # ff = NonlinearFunction(sigma_loop!; sparsity = jac_sparsity)
      # prob = NonlinearSolve.NonlinearProblem(ff, ss_buffer, 0.0)
      prob = NonlinearSolve.NonlinearProblem(sigma_loop!, ss_buffer[left_border:right_border], 0.0)
      sol = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), reltol=1e-6, maxiters=1000)

      # prob = NonlinearSolve.NonlinearLeastSquaresProblem(sigma_loop!, ss_buffer[left_border:right_border], 0.0)
      # sol = NonlinearSolve.solve(prob, Tsit5(), abstol=1e-12, maxiters=1000)
      @. ss_buffer[left_border:right_border] = sol.u
      ######################### END METHOD


      # ## filtering 
      # sigma_freq = fft(ss_buffer)
      # filter::Array{Float64} = ones(M)
      # window = 15 # over 256
      # for i in window:M-window
      #   filter[i] = 0.0
      # end
      # ss_buffer .= real(ifft(sigma_freq .* filter))

      if !all(ss_buffer .>= 0)
        @warn "NEGATIVE sigma values "
        throw(DomainError(-999))
      end
      if !all(ss_buffer .<= 1.01)
        info && @warn "sigma > 1.0"
      end
      for i in 1:M-1
        if (ss_buffer[i+1] - ss_buffer[i]) > 0.05
          info && @warn "discontinuity found"
        end
      end
      # clamp!(ss_buffer, 0.0, 1.0)
      # save solution for next iteration      
      # debug info
      # auxiliary[] = maximum(sol.resid)

      # @info @sprintf("Linf residue= %2.1e" , aux)
      # display(sol.stats)
      # info && @info sol.retcode

    catch err
      if isa(err, DomainError)
        tmp_real1 = NaN
        throw(NpseCollapse(-999))
      else
        throw(err)
      end
    end

    ## tmp_real1 <--- sigma^2 from the solution
    ## tmp_real2 <--- derivative of sigma
    tmp_real1 .= ss_buffer .^ 2
    tmp_real2[1] = (ss_buffer[2] - 1.0) / dxx
    tmp_real2[M] = (1.0 - ss_buffer[M-1]) / dxx
    @inbounds for i = 2:M-1
      tmp_real2[i] = (ss_buffer[i+1] - ss_buffer[i-1]) / dxx
    end
    nonlinear =
      g * abs2.(psi) ./ tmp_real1 + (1 / 2 * tmp_real1) .+
      (1 ./ (2 * tmp_real1)) .* (1 .+ (tmp_real2 .^ 2))
    @. psi = exp(dt_order * -im * iswitch * (V0 + nonlinear)) * psi
  end

  if equation != GPE_1D
    tmp_real1 .= abs2.(psi)
    if maximum(tmp_real1) > collapse_threshold / dV
      throw(Gpe3DCollapse(maximum(abs2.(psi) * dV)))
    end
  end
  kspace!(psi, sim)
  return nothing
end


@inline function propagate_manual!(
  psi::Array{ComplexF64},
  psi_i::Array{ComplexF64},
  tmp_real1::Array{Float64},
  tmp_real2::Array{Float64},
  sim::Sim{1,Array{ComplexF64}},
  t,
  aux;
  ss_buffer::Array{Float64}=nothing,
  info=false,
)
  (ksquared, iswitch, mu, gamma_damp, dt) =
    unpack_selection(sim, :ksquared, :iswitch, :mu, :gamma_damp, :dt)
  # splitting: N/2, N/2, L
  @. psi =
    exp(dt * iswitch * (1.0 - im * gamma_damp) * (-im * (1 / 2 * ksquared - mu))) * psi
  nlin_manual!(psi, tmp_real1, tmp_real2, sim, t, aux; ss_buffer=ss_buffer, info=info)
  # nlin_manual!(psi, tmp_psi2, real_psi, sim, t, aux; ss_buffer=ss_buffer, info=info)
  if iswitch == -im
    psi .= psi / sqrt(nsk(psi, sim))
    cp_diff =
      (chempotk_simple(psi, sim) - chempotk_simple(psi_i, sim)) /
      chempotk_simple(psi_i, sim) / dt
    psi_i .= psi
    return cp_diff
  else
    return 0.0
  end
end

