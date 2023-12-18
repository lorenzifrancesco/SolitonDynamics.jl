# General solvers to be used with DiffEq

function nlin!(dpsi, psi, sim::Sim{1,Array{ComplexF64}}, t)
    @unpack ksquared, g, X, V0, iswitch, dV, Vol, mu, equation, sigma2, N = sim
    x = X[1] |> real
    N = N[1]
    dpsi .= psi
    mu_im = 0.0
    if iswitch == -im
        mu_im = chempotk(psi, sim) # wainting for a more efficient implementation
    end
    xspace!(dpsi, sim)
    if equation == GPE_1D
        @. dpsi *= -im * iswitch * (V0 + V(x, t) + g * abs2(dpsi)) + mu_im
    elseif equation == NPSE
        nonlinear =
            g * abs2.(dpsi) ./ sigma2.(dpsi) +
            (1 ./ (2 * sigma2.(dpsi)) + 1 / 2 * sigma2.(dpsi))
        @. dpsi *= -im * iswitch * (V0 + V(x, t) + nonlinear) + mu_im
    elseif equation == NPSE_plus
        @assert false # still have to implement the corrections
        sigma2_plus = zeros(length(x))
        try
            # Nonlinear Finite Element routine
            b = (1 .+ g * abs2.(dpsi))
            b[1] += 1.0 * 1 / (4 * dV)
            b[end] += 1.0 * 1 / (4 * dV)
            a = ones(length(b))
            A0 = 1 / (2 * dV) * SymTridiagonal(2 * a, -a)
            ss = ones(N)
            prob = NonlinearProblem(sigma_eq, ss, [b, A0, dV])
            sol = solve(prob, NewtonRaphson(), reltol = 1e-3)
            sigma2_plus = (sol.u) .^ 2
        catch err
            if isa(err, DomainError)
                sigma2_plus = NaN
                throw(NpseCollapse(-666))
            else
                throw(err)
            end
        end
        nonlinear =
            g * abs2.(dpsi) ./ sigma2_plus + (
                1 / 2 * sigma2_plus .+
                (1 ./ (2 * sigma2_plus)) .*
                (1 .+ (1 / dV * diff(prepend!(sigma2_plus, 1.0))) .^ 2)
            )
        # warning: sigma2_plus gets modified by prepend!
        @. dpsi *= -im * iswitch * (V0 + V(x, t) + nonlinear) + mu_im
    elseif equation == CQGPE
        throw("unimplemented")
        @. dpsi =
            (-im * iswitch * (V0 + V(x, t) + g * abs2(dpsi))) * dpsi +
            abs2(abs2(dpsi)) * dpsi
    end
    kspace!(dpsi, sim)
    return nothing
end


function propagate!(dpsi, psi, sim::Sim{1,Array{ComplexF64}}, t; info = false)
    @unpack ksquared, iswitch, dV, Vol, mu, gamma_damp = sim
    nlin!(dpsi, psi, sim, t)
    @. dpsi = (1.0 - im * gamma_damp) * (-im * (1 / 2 * ksquared - mu) * psi + dpsi)
    return nothing
end
