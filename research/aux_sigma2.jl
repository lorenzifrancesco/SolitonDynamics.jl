function gs_sigma2()
  gamma_list = [0.65]
  for gamma in gamma_list
    sd = load_parameters_collapse(gamma_param=gamma)
    prepare_for_collision!(sd)

    p=plot(title="gamma=$(gamma)")
    for (k, v) in sd
      est = estimate_sigma2(v.psi, v)
      plot!(p, real(v.X[1]), est)
    end
    display(p)
  end
end
