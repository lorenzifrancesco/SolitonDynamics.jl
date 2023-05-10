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
    axial_density = sum(aa, dims=(2, 3))[:, 1, 1] * sim.dV
    for x in xax
        for y in yax
            for z in zax
                if axial_density[x] < 1e-50
                    tmp[x, y, z] = aa[x, y, z] / axial_density[x] * sim.dV
                    @warn "found small prob"
                else
                    @time tmp[x, y, z] = sim.X[2][y]^2 * sim.X[3][z]^2 * aa[x, y, z] / axial_density[x] * sim.dV
                end
            end
        end
    end

    s2 = 2 * sum(tmp, dims=(2, 3))[:, 1, 1]
    return s2
end

function project_radial(psi_k,sim::Sim{3, CuArray{ComplexF64}})
    # MSE estimator
    psi = xspace(psi_k, sim)
    aa = Array(abs2.(psi))
    x = sim.X[1]
    y = sim.X[2]
    z = sim.X[3]
    rmax = sim.X[2] * sqrt(2)
    r_steps = 256
    r_axis = LinRange(0, rmax, r_steps) |> collect
    dr = rmax / r_steps
    radial_density = zeros((sim.N[1], r_steps)) # axis and radius
    for (ix, x) in eunerate(x)
        for (iy, y) in enumerate(y)
            for (iz, z) in enumerate(z)
                distance = sqrt(y^2 + z^2) 
                ir = floor(r_steps * distance/rmax)
                radial_density[ix, ir] += aa[ix, iy, iz] * dV
            end
        end
    end
    radial_density .= radial_density / dr
    return r_axis, radial_density
end