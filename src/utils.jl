"""
estimate the Gaussian sigma2 parameter from 3D data
"""
function estimate_sigma2(psi,sim::Sim{3, CuArray{ComplexF64}})
    s2 = Array{Float64, 1}(undef, sim.N[1])
    # MSE estimator
    dx = sim.X[1][2] - sim.X[1][1]
    aa = Array(abs2.(psi))
    xax = 1:sim.N[1]
    yax = 1:sim.N[2]
    zax = 1:sim.N[3]
    tmp = zeros(sim.N)
    axial_density = sum(aa, dims=(2, 3))[:, 1, 1] * sim.dV
    @info sum(axial_density * dx)
    @info size(axial_density)
    for x in xax
        for y in yax
            for z in zax
                tmp[x, y, z] = sim.X[2][y]^2 * sim.X[3][z]^2 * aa[x, y, z] / axial_density[x] * sim.dV
            end
        end
    end
    @info size(tmp)
    s2 = sum(tmp, dims=(2, 3))[:, 1, 1]
    return s2
end