function load_parameters_gs(; gamma_param::Float64=0.6, eqs=["G1", "N", "Np", "G3"])
    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    N_axial_steps = 512
    abstol_all = 1e-12
    initial_width = 100 
    dt_all = 0.1
    # =========================================================
    ## 1D-GPE 
    L = (40.0,)
    N = (N_axial_steps,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
    @unpack_Sim sim_gpe_1d

    iswitch = -im
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    maxiters = maxiters_1d
    g = - 2 * gamma_param
    abstol = abstol_all
    dt = dt_all
    x = X[1]
    k = K[1]
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    tf = 1e10
    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param))
    @. psi_0 = exp(-x^2/initial_width)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE plus (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    N_axial_steps = 512
    L_axial = 40.0
    L = (L_axial,10.0,10.0)
    N = (N_axial_steps, 64, 64)
    dx = L_axial / N_axial_steps
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = -im
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4*pi)
    abstol = abstol_all
    maxiters = maxiters_3d
    dt = dt_all / 10
    x0 = 0.0
    vv = 0.0
    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    tf = 1e10
    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*x*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    initial_3d = copy(psi_0)
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end

function load_parameters(; vv::Float64 = 0.0, bb::Float64 = 0.0, gamma_param::Float64=0.6, Nsaves::Int64=200, eqs=["G1", "N", "Np", "G3"])
    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    dt_all = 0.1
    iswitch_all = -im

    max_vel = 1.0
    N_axial_steps = 2048
    abstol_all = 1e-8
    initial_width = 100
    Lx = 40.0
    # =========================================================
    ## 1D-GPE 
    L = (Lx,)
    N = (N_axial_steps,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
    @unpack_Sim sim_gpe_1d

    iswitch = iswitch_all
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    g = - 2 * gamma_param
    n = 100
    abstol = abstol_all
    x = X[1]
    k = K[1]
    x0 = L[1] / 4
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    alg = BS3()

    # will be overwritten
    time_steps = 500
    Nt = Nsaves
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end 
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_1d
    dt = dt_all

    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param)) 
    @. psi_0 = exp(-(x-x0)^2/initial_width) * exp(-im*(x-x0)*vv)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    Nx = 512
    L = (Lx,10.0,10.0)
    N = (Nx, 64, 64)
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = iswitch_all
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4*pi)
    abstol = abstol_all
    alg = BS3()
    
    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    time_steps = 500
    Nt = Nsaves
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    if vv > max_vel/2
        time_steps = 1000
    elseif vv > max_vel/4
        time_steps = 2500
    else
        time_steps = 4000
    end
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_3d
    dt = dt_all

    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*(x-x0)*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end

function load_parameters_alt(; vv::Float64 = 0.0, bb::Float64 = 0.0, gamma_param::Float64=0.6, Nsaves::Int64=200, eqs=["G1", "N", "Np", "G3"])
    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    dt_all = 0.01
    iswitch_all = -im

    max_vel = 1.0
    N_axial_steps = 1024
    abstol_all = 1e-4
    initial_width = 10
    Lx = 80.0
    # =========================================================
    ## 1D-GPE 
    L = (Lx,)
    N = (N_axial_steps,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
    @unpack_Sim sim_gpe_1d

    iswitch = iswitch_all
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    g = - 2 * gamma_param
    n = 100
    abstol = abstol_all
    x = X[1]
    k = K[1]
    x0 = L[1] / 4
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    alg = BS3()

    # will be overwritten
    time_steps = 500
    Nt = Nsaves
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end 
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_1d
    dt = dt_all

    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param)) 
    @. psi_0 = exp(-(x-x0)^2/initial_width) * exp(-im*(x-x0)*vv)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    Nx = 512
    L = (Lx,10.0,10.0)
    N = (Nx, 64, 64)
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = iswitch_all
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4 * pi)
    abstol = abstol_all
    alg = BS3()
    # we can augment the accuracy


    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    time_steps = 500
    Nt = Nsaves
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    if vv > max_vel/2
        time_steps = 1000
    elseif vv > max_vel/4
        time_steps = 2500
    else
        time_steps = 4000
    end
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_3d
    dt = dt_all

    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*(x-x0)*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end

function prepare_in_ground_state!(sim::Sim{1, Array{Complex{Float64}}})
    # compute the ground state
    # start from a convenient initial state (it doesn't matter by the way)
    @unpack_Sim sim
    x0 = L[1] / 4
    iswitch = -im 
    maxiters = 1e10
    abstol = 1e-6
    initial_width = 5
    @. psi_0 = exp(-((X[1]-x0)/initial_width)^2/2)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    x = X[1] |> real
    V0 = zeros(N[1])
    tmp_dt = deepcopy(dt)
    dt = 0.005
    @pack_Sim! sim

    @info "Computing ground state..."
    sol = runsim(sim; info=false)

    @info "Assigning GS as dynamical sim initial state..."
    xspace!(sol.u, sim)
    display(sol.u)

    # pack everything back up, imprint the correct velocity (suppose x0 stays constant)
    @unpack_Sim sim
    iswitch = 1
    x = X[1]
    dt = tmp_dt
    @. psi_0 = sqrt(abs2(sol.u))
    kspace!(psi_0, sim)
    @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
    @pack_Sim! sim
    return sim
end

function prepare_in_ground_state!(sim::Sim{3, CuArray{Complex{Float64}}})
    # compute the ground state
    # start from a convenient initial state (it doesn't matter by the way)
    @unpack_Sim sim
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    x0 = L[1] / 4
    iswitch = -im 
    maxiters = 1e10
    abstol = 1e-8
    initial_width = 5
    tmp_dt = deepcopy(dt)
    dt = 0.005
    tmp = [exp(-(((x-x0)/initial_width)^2 /2 + (y^2 + z^2)/2)) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    V0 *= 0.0
    psi_0 = kspace(psi_0, sim) # FIXME strangest behaviour ???!
    @pack_Sim! sim

    @info "Computing ground state..."
    sol = runsim(sim; info=true)

    @info "Assigning GS as dynamical sim initial state..."
    xspace!(sol.u, sim)
    print("size of sol.u: ", size(sol.u))

    # pack everything back up, imprint the correct velocity (suppose x0 stays constant)
    @unpack_Sim sim
    iswitch = 1
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    dt = tmp_dt
    @. psi_0 = sqrt(abs2.(sol.u))
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
    @pack_Sim! sim
    return sim
end

function imprint_vel_set_bar(sim::Sim{1, Array{Complex{Float64}}}; vv::Float64=0.0, bb::Float64=0.0, bw::Float64=0.5)
    simc = deepcopy(sim)
    @unpack_Sim simc
    x = X[1] |> real
    @. V0 = bb * exp(-(x/bw)^2 /2) # central barrier
    x0 = L[1]/4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps
    xspace!(psi_0, simc)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, simc)
    @pack_Sim! simc
    return simc
end

function imprint_vel_set_bar(sim::Sim{3, CuArray{Complex{Float64}}}; vv::Float64=0.0, bb::Float64=0.0, bw::Float64=0.5)
    simc = deepcopy(sim)
    @unpack_Sim simc
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    V0 = [1/2*(z^2 + y^2) + bb * exp(-(x/bw)^2 /2) for x in x, y in y, z in z] # central barrier
    x0 = L[1]/4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps
    xspace!(psi_0, simc)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, simc)
    @pack_Sim! simc
    return simc
end

function imprint_vel_set_bar!(sim::Sim{1, Array{Complex{Float64}}}; vv::Float64=0.0, bb::Float64=0.0, bw::Float64=0.5)
    @unpack_Sim sim
    x = X[1] |> real
    @. V0 = bb * exp(-(x/bw)^2 /2) # central barrier
    x0 = L[1]/4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps
    xspace!(psi_0, sim)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, sim)
    @pack_Sim! sim
    return sim
end

function imprint_vel_set_bar!(sim::Sim{3, CuArray{Complex{Float64}}}; vv::Float64=0.0, bb::Float64=0.0, bw::Float64=0.5)
    @unpack_Sim sim
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    x0 = L[1]/4
    V0 = [1/2*(z^2 + y^2) + bb * exp(-(x/bw)^2 /2) for x in x, y in y, z in z] # central barrier
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps
    xspace!(psi_0, sim)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, sim)
    @pack_Sim! sim
    return sim
end

function set_g!(sim::Sim{1, Array{Complex{Float64}}}, gamma_param::Float64=0.4)
    @unpack_Sim sim
    g = -2 * gamma_param
    @pack_Sim! sim
    return
end

function set_g!(sim::Sim{3, CuArray{Complex{Float64}}}, gamma_param::Float64=0.4)
    @unpack_Sim sim
    g = - (4 * pi) * gamma_param
    @pack_Sim! sim
    return
end