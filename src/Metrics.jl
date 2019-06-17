module Metrics

import DelimitedFiles: writedlm, readdlm

import CellSimCommon
const CSC = CellSimCommon

import Cortex
import Nucleus

import YAML

function init_metrics(iter::Int64,
                      P::CSC.Params, F::CSC.Flags, config,
                      coords::Cortex.PointCoords,
                      nucleus_coords::Union{Nucleus.NucleusCoords,Missing})
    m = Dict{String, Union{Int64, Float64, Vector{Float64}}}()
    m["iter_start"] = iter
    m["t_start"] = P.δt * iter
    m["max_y_start"] = maximum(coords.x[:,2])
    m["target_max_y"] = m["max_y_start"] + config["metrics"]["periods"]*2π/P.f_ω0

    m["inst_max_y"] = Vector{Float64}()
    m["inst_velocity"] = Vector{Float64}()
    m["inst_cortex_area"] = Vector{Float64}()
    m["inst_cortex_perimeter"] = Vector{Float64}()
    if F.nucleus
        m["inst_nucleus_area"] = Vector{Float64}()
        m["inst_nucleus_perimeter"] = Vector{Float64}()
        if F.centrosome
            m["inst_n2c_distance"] = Vector{Float64}()
        end
    end
    return m
end

function update_metrics(m::Dict, iter::Int64,
                        P::CSC.Params, F::CSC.Flags, config,
                        coords::Cortex.PointCoords,
                        coords_s::Cortex.PointCoordsShifted,
                        coords_old::Cortex.PointCoords,
                        nucleus_coords::Union{Missing, Nucleus.NucleusCoords})

    push!(m["inst_max_y"], maximum(coords.x[:,2]))
    push!(m["inst_velocity"], (m["inst_max_y"][end] - maximum(coords_old.x[:,2]))/P.δt)

    area = 0.25*sum(coords.ΔL .* CSC.@dotprod(coords.v, coords.x) + coords_s.ΔL_p .* CSC.@dotprod(coords_s.v_p, coords.x))
    push!(m["inst_cortex_area"], area)
    push!(m["inst_cortex_perimeter"], sum(coords.ΔL))

    if F.nucleus
        c = nucleus_coords
        nucleus_area = 0.25*sum(c.q.*(CSC.@dotprod(c.n, c.Y) + CSC.@dotprod(c.n[c.circ_idx.p1,:], c.Y)))
        push!(m["inst_nucleus_area"], nucleus_area)
        push!(m["inst_nucleus_perimeter"], sum(c.r))
        if F.centrosome
            push!(m["inst_n2c_distance"], sqrt(sum(abs2, coords.centro_x - vec(sum(nucleus_coords.Y; dims=1))/P.Nnuc)))
        end
    end
end

function close_metrics(m::Dict, iter::Int64, P::CSC.Params)
    m["t_end"] = iter*P.δt
    m["max_y_end"] = m["inst_max_y"][end]
    m["average_velocity"] = (m["max_y_start"] - m["max_y_end"])/(m["t_end"] - m["t_start"])
end

function save_metrics(m::Dict, prefix::String)
    f = open(string(prefix, "_metrics.yaml"), "w")
    for k in ["iter_start", "max_y_start", "t_start", "target_max_y", "max_y_end", "t_end", "average_velocity"]
        write(f, string(k, ": ", m[k], "\n"))
    end
    close(f)

    for k in ["inst_max_y", "inst_velocity", "inst_cortex_area", "inst_cortex_perimeter", "inst_nucleus_area", "inst_nucleus_perimeter", "inst_n2c_distance"]
        writedlm(string(prefix, "_metrics_", k, ".csv"), m[k], ',')
    end
end

function read_metrics(prefix::String)
    m = YAML.load(open(string(prefix, "_metrics.yaml"), "r"))

    for k in ["inst_max_y", "inst_velocity", "inst_cortex_area", "inst_cortex_perimeter", "inst_nucleus_area", "inst_nucleus_perimeter", "inst_n2c_distance"]
        m[k] = readdlm(string(prefix, "_metrics_", k, ".csv"), ',')
    end

    return m
end


end # module
