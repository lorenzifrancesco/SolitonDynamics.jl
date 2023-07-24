function compare_chempot(; use_precomputed=false)
  @info "Loading parameters..."
  sd = load_parameters_alt(eqs=["G1", "N"])
  FFTW.set_num_threads(1)
  gamma_range = LinRange(0.65, 0.75, 30)
  p = plot()

  if isfile("results/mu_db.jld2")
    mud = load("results/mu_db.jld2")
  else
    mud = Dict()
    save("results/mu_db.jld2", "mud", mud)
  end
  for (k, sim) in sd
    @info "==>Running for $k"

    if haskey(mud, hs(k, 0.666)) && use_precomputed
      @info "--> Loading for $k..."
      sol = mud[k]
    else
      @info "--> Computing for $k..."

      mu = zeros(length(gamma_range))
      for (ig, gamma) in enumerate(gamma_range)
        set_g!(sim, gamma)
        print("\n>>>Computing for gamma = $gamma")
        sol = nothing
        try
          sol = runsim(sim; info=false)
        catch e
          if isa(e, NpseCollapse) || isa(e, Gpe3DCollapse)
            @warn "Collapse detected for gamma = $gamma"
            mu[ig] = NaN
            continue
          else
            rethrow(e)
          end
        end
        print("-->final time = $(sol.cnt * sim.dt) \n --> Computing chempot...")
        mu[ig] = chempotk(sol.u, sim)
        print("mu = $(mu[ig])\n")
      end

      push!(mud, hs(k, 0.666) => mu)
      save("result/mu_db.jld2", "mud", mud)
    end
    plot!(p, gamma_range, mu, label=nameof(k), color=colorof(k), linestyle=lineof(k))
  end
  savefig(p, "media/chempot_compare.pdf")
end