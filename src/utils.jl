"""
estimate the Gaussian sigma2 parameter from 3D data
"""
function estimate_sigma2(psi,sim::Sim{3, CuArray{ComplexF64}})
    s2 = Array{Float64, 1}(undef, sim.N[1])
    # MSE estimator
    distance2 = Array(abs2.(sim.X[2])) * Array(ones(sim.N[2]))' .+ Array(abs2.(sim.X[3])) * Array(ones(sim.N[3]))'
    aa = Array(abs2.(psi))
    xax = 1:sim.N[1]
    yax = 1:sim.N[2]
    zax = 1:sim.N[3]
    @info xax
    tmp = [distance2[y, z] * aa[x, y, z] for x in xax for y in yax for z in zax]
    s2 = sum(tmp, dims=(2, 3))[:, 1, 1]
    @info size(s2)

    return s2
end