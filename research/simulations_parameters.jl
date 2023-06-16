using CondensateDynamics
using OrdinaryDiffEq
import FFTW
using CUDA
# import HDF5

function get_parameters()
    maxiters_1d = 1e10
    maxiters_3d = 1e10
    N_axial_steps = 1024
    abstol_all = 1e-8
    gamma_param = 0.15 
    initial_width = 100 
    
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
    g_param = gamma_param
    maxiters = maxiters_1d
    g = - 2 * g_param
    n = 100
    as = g_param / n
    abstol = abstol_all
    dt = 0.005
    x = X[1]
    k = K[1]
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    width = 7
    tf = 1e10
    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(g_param/2) * 2/(exp(g_param*x) + exp(-x*g_param))
    psi_0 .= exp.(-x.^2/initial_width)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    # =========================================================
    ## NPSE (unable to copy)

    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if g_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if g_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus
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
    g = - g_param * (4*pi)
    abstol = abstol_all
    maxiters = maxiters_3d
    dt = 0.01
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
    # =========================================================
    sim_dictionary = Dict(
                        "GPE_1D_GS" => sim_gpe_1d, 
                        "NPSE_GS" => sim_npse, 
                        "NPSE_plus_GS" => sim_npse_plus, 
                        "GPE_3D_GS" => sim_gpe_3d,
                        )

    return sim_dictionary
end