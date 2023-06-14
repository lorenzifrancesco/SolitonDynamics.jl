"""
estimate the Gaussian sigma2 parameter from 3D data
"""
function estimate_sigma2(psi_k,sim::Sim{3, CuArray{ComplexF64}})
    s2 = Array{Float64, 1}(undef, sim.N[1])
    # MSE estimator
    psi = xspace(psi_k, sim)
    dx = sim.X[1][2] - sim.X[1][1]
    aa = Array(abs2.(psi))
    xax = 1:sim.N[1]
    yax = 1:sim.N[2]
    zax = 1:sim.N[3]
    tmp = zeros(sim.N)
    axial_density = sum(aa, dims=(2, 3))[:, 1, 1]
    for x in xax
        for y in yax
            for z in zax
                if axial_density[x] < 1e-50
                    tmp[x, y, z] = aa[x, y, z] / axial_density[x]
                    @warn "found small prob"
                else
                    tmp[x, y, z] = sim.X[2][y]^2 * sim.X[3][z]^2 * aa[x, y, z] / axial_density[x]
                end
            end
        end
    end
    s2 = 4 * sum(tmp, dims=(2, 3))[:, 1, 1]
    return s2
end

function project_radial(psi_k,sim::Sim{3, CuArray{ComplexF64}})
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
                ir = Int(round(r_steps * distance/rmax))
                radial_density[ix, ir] += aa[ix, iy, iz] * sim.dV
            end
        end
    end
    #radial_density[:, :] = aa[:, :, 64]
    radial_density .= radial_density / dr
    return r_axis, radial_density
end


function npse_expr(mu) 
    f =  ((1-mu)^(3/2) - 3/2*(1-mu)^(1/2))*(2*sqrt(2))/3
    return f
end

"""
Compute the ground state energy normalized to the harmonic energy unit
"""
function npse_energy(n, as)
    steps = 500
    dn = n/steps
    energy = 0
    for n_int in LinRange(0, n, steps)
        gamma = as * n_int
        a = roots(x->npse_expr(x)+gamma, 0.5..1, Newton, 1e-4)
        mu = mid(a[1].interval)
        energy += mu
    end
    return energy * dn
end

function gpe_energy(n, as)
    steps = 500
    dn = n/steps
    energy = 0
    for n_int in LinRange(0, n, steps)
        gamma = as * n_int
        mu = 1-gamma^2/2
        energy += mu
    end
    return energy * dn
end

"""
Compute the chemical potential
"""
function npse_mu(n, as)
    gamma = as * n
    a = roots(x->npse_expr(x)+gamma, 0.5..1, Newton, 1e-4)
    mu = mid(a[1].interval)
    return mu - 1
end

function gpe_mu(n, as)
    return - (n * as)^2/8
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
    dV,dk = measures(sim.L, sim.N)
    return sum(abs2.(psi)) * dk
end

"""
compute norm squared of a region
"""
function ns(psi, sim, mask)
    return sum(abs2.(psi).* mask) * sim.dV
end

"""
return σ^2(ψ) of the NPSE
"""
function init_sigma2(g::Float64)
    function sigma2(psi::ComplexF64)
        result = psi
        try
            result = sqrt(1 + g * abs2(psi))
        catch  err
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

# """
# check notes, Sigma is sigma^2
# """
# function deSigmaLambda!(dsigma, sigma, psi, t)
#     return 
# end

# """
# return the sigma2 function to be called in the ODE loop
# """
# function init_sigma2_ode(g::Float64, x::Array{Float64})
#     function sigma2(psi::Array{ComplexF64}, sigma2_0::Float64, lambda::Float64)
#         problem = ODEProblem(deSigmaLambda!, Array([sigma2_0, lambda]), (x[1], x[end]), psi)
#         sol = solve(problem,
#                     alg=BS3(),
#                     dense=false,
#                     maxiters=2000,
#                     progress=true, 
#                     )
#         result = sol.u
#         try
#             tmp = sqrt.(result)
#         catch  err
#             if isa(err, DomainError)
#                 result = NaN
#                 throw(NpseCollapse(NaN))
#             else
#                 throw(err)
#             end
#         end
#         return result
#     end
#     return sigma2
# end

"""
chemical potential in a given configuration
"""
function chempotk(psi, sim)
    @unpack ksquared,dV,V0,Vol,g = sim 
    mu = 1/Vol * sum(1/2 * ksquared .* abs2.(psi))
    tmp = xspace(psi, sim)
    mu += dV * sum((V0 + g*abs2.(tmp)) .* abs2.(tmp))
    mu *= 1/ns(tmp, sim) 
    #mu += 1 # add one transverse energy unit (1D-GPE case)
    return mu
end

"""
chemical potential in a given configuration
"""
function chempot(psi, sim)
    @unpack ksquared,dV,V0,Vol,g = sim 
    mu = dV * sum((V0 + g*abs2.(psi)) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1/Vol * sum(1/2 * ksquared .* abs2.(tmp))
    mu *= 1/nsk(tmp, sim) 
    #mu += 1 # add one transverse energy unit (1D-GPE case)
    return mu
end