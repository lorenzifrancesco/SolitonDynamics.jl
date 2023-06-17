using BenchmarkTools
using CondensateDynamics

function run_benchmarks()
    include("research/simulations_parameters.jl")
    simulations_dict = get_parameters()
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