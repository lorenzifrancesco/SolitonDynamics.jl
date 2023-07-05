using ColorSchemes

function gs_sigma2()
  gamma_list = [0.55]
  for gamma in gamma_list
    sd = load_parameters_alt(gamma_param=gamma)
    prepare_for_collision!(sd, gamma; use_precomputed_gs=true, info=true)

    p=plot(title="gamma=$(gamma)")
    if length(sd)>1
      pal = palette([:red, :blue], length(sd))
    else
      pal = [:red]
    end
    i = 1
    for (k, v) in sd
      plot_axial_density(v.psi_0, v; show=true)
      est = estimate_sigma2(v.psi_0, v)
      plot!(p, real(v.X[1]), est, color=pal[i], label=k)
      # mid = Int(round(v.N[1]/2))
      # mid_tran = Int(round(v.N[2]/2))
      # pp = xspace(v.psi_0, v)
      # tran = real(v.X[2])
      # dx = real(v.X[1][2] - v.X[1][1])
      # dy = real(v.X[2][2] - v.X[2][1])
      # dz = real(v.X[3][2] - v.X[3][1])
      # @warn sum(abs2.(pp[mid, :, mid_tran])/sum(abs2.(pp[mid, :, mid_tran])))
      # q = plot(tran, abs2.(pp[mid, :, mid_tran])/sum(abs2.(pp[mid, :, mid_tran])) / dx)
      # plot!(q, tran, abs2.(pp[mid, mid_tran, :])/sum(abs2.(pp[mid, mid_tran, :])) / dx)
      # plot!(q, tran, exp.(-tran.^2) / sqrt(pi), color=:red, label="Gaussian")
      # plot!(q, tran, exp.(-tran.^2/1.0212) / sqrt(pi*1.0212), ls=:dot, label="Gaussian with correct width")
      # display(q)
      # show_slice(0.0, v.psi_0, v; show=true)
      # show_profile(0.0, v.psi_0, v; show=true)
      i +=1
    end
    display(p)
  end
end