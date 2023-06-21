includet("CDResearch.jl")


static_standard = load_parameters_gs(gamma_param=0.6)
dynamic_standard = load_parameters_dy(gamma_param=0.6)

if Threads.nthreads() == 1
    @warn "running in single thread mode!"
else
    @info "running in multi-thread mode: n_threads =" Threads.nthreads()
end

# simulations are already loaded 
@info "Static simulations: " keys(static_standard)
@info "Dynamic simulations: " keys(dynamic_standard)

# prepare ground states (saving them)
if isfile("research/sd_syn_prepared.jld2")
    @info "Loading prepared simulations..."
    sd = JLD3.load("research/sd_syn_prepared.jld2", "sd")
else
    @info "Preparing dynamical simulations in the ground state..."
    ksd = keys(dynamic_standard)
    sl = prepare_in_ground_state!.(values(dynamic_standard))
    sd = Dict(zip(ksd, sl))
    JLD3.save("research/sd_syn_prepared.jld2", "sd", sd)
end

# TODO here we can obtain the plots of the GS using the prepared simulations

## ========= run standardized simulations 
## ========= for publication-grade plots

# # -Ground states 
# pgs = all_ground_states()

# Lines
pli = all_lines()

# Tiles (choose a simulation)
tiles = all_tiles()