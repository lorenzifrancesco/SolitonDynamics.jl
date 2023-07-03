using ColorSchemes

function gs_sigma2()
  gamma_list = [0.4]
  for gamma in gamma_list
    sd = load_parameters_collapse(gamma_param=gamma, eqs=["Np"])
    prepare_for_collision!(sd, gamma)

    p=plot(title="gamma=$(gamma)")
    pal = palette([:red, :blue], length(sd))
    i = 1
    for (k, v) in sd
      est = estimate_sigma2(v.psi_0, v)
      plot!(p, real(v.X[1]), est, color=pal[i], label=k)
      i +=1
    end
    display(p)
  end
end