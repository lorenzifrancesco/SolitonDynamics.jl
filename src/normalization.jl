
const hbar = 1.054571817e-34
const amu = 1.6605390666e-27

@with_kw mutable struct SISim{D} <: Simulation{D}
    as::Float64 = -0.2
    omega_perp::Float64 = 1700 * 2 * pi  # Li-7 mass
    mass::Float64 = 7.016 * amu
    N_particles::Int32 = 1e4
    L::NTuple{D,Float64} = (10e-6, 10e-6, 10e-6) # length scales
    N::NTuple{D,Int64} = (256, 256, 256)  # grid points in each dimensions
end

function printsim(sisim::SISim)
    @unpack as, omega_perp, mass, N_particles, L, N = sisim
    @info "Simulation in SI units:"
    Base.print("\nas = ", as)
    sim = normalize(sisim)
    Base.print("\ng = ", sim.g)
end

"""
    normalize(parameters)    
conversion tool from SI units used in experiments to simulation Parameters
"""
function normalize(sisim::SISim)
    @unpack as, omega_perp, mass, N_particles, L, N = sisim
    # choice of units
    time_unit = (2 * pi) / omega_perp # inverse of frequency
    energy_unit = hbar * omega_perp
    length_unit = sqrt(hbar / (mass * omega_perp))
    sim = Sim{3}(
        g = 4 * pi * as * (N_particles - 1) / (mass * hbar * omega_perp),
        L = L ./ length_unit,
        N = N,
    )
    return sim
end
