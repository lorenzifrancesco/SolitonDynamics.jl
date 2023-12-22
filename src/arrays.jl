function xspace(ϕ::AbstractArray, sim::Sim)
    return sim.T.Tkx * ϕ
end

function xspace!(ψ, sim)
    @unpack T = sim
     T.Tkx! * ψ
    return nothing
end

function kspace(ψ, sim)
    @unpack T = sim
    return T.Txk * ψ
end

function kspace!(ψ, sim)
    @unpack T = sim
    T.Txk! * ψ
    return nothing
end
