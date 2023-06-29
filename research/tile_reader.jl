include("plotter_all.jl")

function view_all_tiles()
    tile_file = "results/tile_dict.jld2"
    @assert isfile(tile_file)
    td = load(tile_file)
    for (k, v) in td
        @info "found" ihs(k)
        ht = heatmap(v)
        savefig(ht, "media/" * string(ihs(k)) *  "_tiles.pdf")
        #  display(ht)
    end
end