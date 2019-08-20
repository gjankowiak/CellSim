function save_state(prefix::String, P::CSC.Params, F::CSC.Flags, coords::Cortex.PointCoords, nucleus_coords::Union{Nucleus.NucleusCoords, missing})
    if F.centrosome
        writedlm(string(prefix, "_cortex.csv"), coords.x[1:end-1,:], ',')
        writedlm(string(prefix, "_centrosome.csv"), coords.x[end,:], ',')
    else
        writedlm(string(prefix, "_cortex.csv"), coords.x, ',')
    end

    if F.nucleus
        writedlm(string(prefix, "_nucleus_Y.csv"), nucleus_coords.Y, ',')
        writedlm(string(prefix, "_nucleus_α.csv"), nucleus_coords.α, ',')
        writedlm(string(prefix, "_nucleus_β.csv"), nucleus_coords.β, ',')
        writedlm(string(prefix, "_nucleus_k.csv"), nucleus_coords.k, ',')
        writedlm(string(prefix, "_nucleus_r.csv"), nucleus_coords.r, ',')
        writedlm(string(prefix, "_nucleus_q.csv"), nucleus_coords.q, ',')
        writedlm(string(prefix, "_nucleus_θ.csv"), nucleus_coords.θ, ',')
        writedlm(string(prefix, "_nucleus_n.csv"), nucleus_coords.n, ',')
        writedlm(string(prefix, "_nucleus_η.csv"), nucleus_coords.η, ',')
    end
end

function load_state(P::CSC.Params, F.CSC.Flags, prefix::String)
    x = readdlm(string(prefix, "_cortex.csv"), ',')

    if F.centrosome
        coords, coords_s = Cortex.new_PointCoords(x[1:end-1,:], P)
        coords.centro_x[:] = x[end,:]
    else
        coords, coords_s = Cortex.new_PointCoords(x, P)
    end

    if F.nucleus
        Nucleus.NucleusCoords(
        readdlm(string(prefix, "_nucleus_Y.csv"), ','),
        readdlm(string(prefix, "_nucleus_α.csv"), ','),
        readdlm(string(prefix, "_nucleus_β.csv"), ','),
        readdlm(string(prefix, "_nucleus_k.csv"), ','),
        readdlm(string(prefix, "_nucleus_r.csv"), ','),
        readdlm(string(prefix, "_nucleus_q.csv"), ','),
        readdlm(string(prefix, "_nucleus_θ.csv"), ','),
        readdlm(string(prefix, "_nucleus_n.csv"), ','),
        readdlm(string(prefix, "_nucleus_η.csv"), ','),
        0.0,
        CSC.init_circ_idx(P.Nnuc)
        )
    else
        return (coords, missing)
    end
end
