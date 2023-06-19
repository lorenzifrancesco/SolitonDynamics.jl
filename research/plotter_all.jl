includet("CDResearch.jl")

if Threads.nthreads() == 1
    @warn "running in single thread mode!"
else 
    @info "running in multi-thread mode: n_threads =" Threads.nthreads()
end

# simulations are already loaded 
@info "Static simulations: " keys(static_standard)
@info "Dynamic simulations: " keys(dynamic_standard)

## ========= run standardized simulations 
## ========= for publication-grade plots

# Ground states 
pgs = all_ground_states()

# Lines
pli = all_lines()

# Tiles (choose a simulation)
tiles = all_tiles()