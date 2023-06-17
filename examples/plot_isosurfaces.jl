# ================ plotting functions
function dense(phi)
    psi = xspace(phi,sim)
    density = abs2.(psi)
    pmax = maximum(density)
    if pmax == 0
        throw("Maximum density is null")
    end
    return density/pmax
end

function isosurface_animation(sol,Nt, sim;file="3Devolution.gif",framerate=3)
    saveto=joinpath("media",file)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    iter = [Array(xspace(sol[k], sim)) for k in 1:Nt]
    iter = [abs2.(iter[k]) for k in 1:Nt]
    fig = Makie.volume(Makie.@lift(iter[$tindex]/maximum(iter[$tindex])),
                        algorithm =:iso,
                        isovalue=0.2,
                        isorange=0.1,
                        transparency=true
    )

    R = 180
    Makie.record(fig, saveto, 1:Nt; framerate=framerate) do i
        tindex[] = i
    end
    return
end

function isosurface(sol)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    psol = Array(abs2.(xspace(sol, sim)))
    scene = Makie.volume(psol/maximum(psol),
                        algorithm =:iso,
                        isovalue=0.1,
                        isorange=0.1,
    )
    display(scene)
    return
end