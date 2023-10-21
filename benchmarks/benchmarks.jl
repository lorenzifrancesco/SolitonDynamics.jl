using BenchmarkTools
using CondensateDynamics

function run_benchmarks()
    include("research/init.jl")
    simulations_dict = load_parameters_gs()
    @info "Available simulations" simulations_dict
    bms = BenchmarkTools.Trial[]
    for (sim_name, sim) in simulations_dict
        @info "Benchmarking $sim_name"
        bb = @benchmark runsim($sim)
        display(bb)
        push!(bms, bb)
    end
    # WISH implement also dynamics (most crucial part)
    return bms
end