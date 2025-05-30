__precompile__(false)
module Plotting

import Printf

import PyPlot
import PyCall

animation = PyCall.pyimport("matplotlib.animation")

import CellSimCommon
const CSC = CellSimCommon
import Centrosome
import Cortex
import Nucleus
import Wall
import JankoUtils

function init_plot(coords::Cortex.PointCoords, P::CellSimCommon.Params, F::CellSimCommon.Flags)

    PyPlot.ion()

    x = coords.x

    fig = PyPlot.figure(figsize=(12.8, 10), dpi=50)

    if PyPlot.matplotlib.backends.backend != "agg"
        PyPlot.show()
    end
    # figManager = PyPlot.get_current_fig_manager()
    # figManager.window.showMaximized()

    ax = PyPlot.axes(xlim = (-5, 5),ylim=(-15, 15))

    ax.set_aspect("equal", "datalim")

    # Cortex
    ax.plot(x[:,1], x[:,2], ".-", zorder=20)[1]

    # Cortex fillin
    ax.fill(x[:,1], x[:,2], color="#f713e033", zorder=10)[1]

    if F.nucleus
        # Nucleus
        ax.plot(zeros(P.Nnuc), zeros(P.Nnuc), ".-", color="yellow", zorder=16)

        # Nucleus fillin
        ax.fill(zeros(P.Nnuc), zeros(P.Nnuc), color="xkcd:light yellow", zorder=16)[1]

        # Normals
        # ax.quiver(zeros(P.Nnuc), zeros(P.Nnuc),
                    # zeros(P.Nnuc), zeros(P.Nnuc), zorder=100, units="xy", scale_units="xy", scale=1e1, width=1e-2)
        # Velocity
        # ax.quiver(zeros(P.Nnuc), zeros(P.Nnuc),
                    # zeros(P.Nnuc), zeros(P.Nnuc), color="blue", zorder=100, units="xy", scale_units="xy", scale=1e1, width=5e-3)
    end

    # Drag force
    if F.plot_drag
        ax.scatter(x[:,1], x[:,2], color="black", zorder=30)
    end

    # Initial condition
    # ax.plot(x[:,1], x[:,2], color="black", lw=0.5, zorder=1)[1]

    # Transport force
    # ax.quiver(coords.x[1], coords.x[2],
                # zeros(size(x,1)), zeros(size(x,1)), zorder=100, units="xy", scale=30, width=0.001, color="blue")

    if F.centrosome
        # Centrosome
        ## Orientation
        ax.quiver(coords.centro_x[1], coords.centro_x[2],
                    [cos.(coords.centro_angle)], [sin.(coords.centro_angle)], zorder=100, units="xy", scale=10, width=0.01, color="#12b600")

        ## Visibility
        ax.fill(x[:,1], x[:,2], color="#32d600", zorder=15, alpha=1.0)[1]

        ## Center
        mt_circle = PyPlot.matplotlib.patches.Circle((coords.centro_x[1], coords.centro_x[2]), Centrosome.MT_RADIUS, zorder=110, alpha=0.5)
        ax.add_artist(mt_circle)

        # Microtubule force
        # ax.quiver(coords.centro_x[1], coords.centro_x[2],
                    # [0], [0], zorder=100, units="xy", scale=1e-2, width=0.01) # MT force

        # ax.quiver([0.0], [0.0], [0], [0], zorder=100, units="xy", scale=5, width=0.001) # MT force

    end

    if F.circular_wall
        y = collect(range(-1; stop=1, length=1000))
        wall_1 = Wall.compute_walls(y, P, F)
        wall_2 = Wall.compute_walls(y, P, F; right=true)
        ax.plot(wall_1[:,1], wall_1[:,2], wall_2[:,1], wall_2[:,2], color="black", lw=0.5)

        levelset_1 = Wall.compute_walls(y, P, F, 2e-4)
        levelset_2 = Wall.compute_walls(y, P, F, 2e-4; right=true)
        ax.plot(levelset_1[:,1], levelset_1[:,2], levelset_2[:,1], levelset_2[:,2], color="red", lw=0.5)
    else

        y = collect(range(-40; stop=40, length=1000))
        wall = Wall.compute_walls(y, P, F)
        ax.plot(wall[:,1], wall[:,2], -wall[:,1], wall[:,2], color="black", lw=0.5)
        levelset = Wall.compute_walls(y, P, F, 2e-4)
        ax.plot(levelset[:,1], levelset[:,2], -levelset[:,1], levelset[:,2], color="red", lw=0.5)
    end

    # Added mass
    # ax.scatter(x[:,1], x[:,2], color="red")

    ax.axvline(0)

    lims = ax.get_xlim()
    x_min, x_max = minimum(x[:,1]), maximum(x[:,1])
    x_mid = 0.5(x_max+x_min)
    x_span = (x_max-x_min)
    y_min, y_max = minimum(x[:,2]), maximum(x[:,2])
    y_mid = 0.5(y_max+y_min)
    y_span = (y_max-y_min)

    ax.set_xlim((x_mid-x_span, x_mid+x_span))
    ax.set_ylim((y_mid-y_span, y_mid+y_span))

    PyPlot.draw()
    # fig.tight_layout()

    sleep(0.001)
    # println("Finishing plot init...")
    # sleep(1)
    return fig
end

function update_plot(coords::Cortex.PointCoords, nucleus_coords::Union{Nucleus.NucleusCoords,Missing}, iter::Int64, time::Real, P::CellSimCommon.Params, F::CellSimCommon.Flags, initializing::Bool, plotables::CellSimCommon.Plotables,
                     vr::Union{Centrosome.VisibleRegion,Missing})

    x = coords.x

    ax = PyPlot.gca()
    if initializing
        prefix = "[INIT]"
    else
        prefix = ""
    end
    ax.set_title(Printf.@sprintf("%s N: %d, iter: %d, T=%fs", prefix, P.N, iter, time))

    lines = ax.lines
    scatters = ax.collections
    patches = ax.patches
    artists = ax.artists

    idx_l = idx_s = idx_p = idx_a = 1

    # Cortex
    lines[idx_l].set_data(x[:,1], x[:,2])
    idx_l += 1

    # Cortex fillin
    patches[idx_p].set_xy(x)
    idx_p += 1

    if F.nucleus
        # Nucleus
        lines[idx_l].set_data([nucleus_coords.Y[:,1]; nucleus_coords.Y[1,1]],
                                [nucleus_coords.Y[:,2]; nucleus_coords.Y[1,2]])
        idx_l += 1

        patches[idx_p].set_xy(nucleus_coords.Y)
        idx_p += 1

        # normal-s
        # midpoints = 0.5*(nucleus_coords.Y+circshift(nucleus_coords.Y, 1))
        # scatters[idx_s].set_offsets(midpoints)
        # scatters[idx_s].set_UVC(nucleus_coords.n[:,1], nucleus_coords.n[:,2])
        # idx_s += 1


        # velocity
        #
        # velocity_x = -nucleus_coords.α. .* nucleus_coords.n[:,2] + nucleus_coords.β. .* nucleus_coords.n[:,1]
        # velocity_y = nucleus_coords.α. .* nucleus_coords.n[:,1] + nucleus_coords.β. .* nucleus_coords.n[:,2]

        # scatters[idx_s].set_offsets(midpoints)
        # scatters[idx_s].set_UVC(velocity_x, velocity_y)
        # idx_s += 1
    end

    # Drag force
    if F.plot_drag
        scatters[idx_s].set_sizes(10CellSimCommon.@entry_norm(plotables.drag_force))
        scatters[idx_s].set_offsets(x)
        idx_s += 1
    end

    # Initial condition
    # idx_l += 1

    # Transport force
    # scatters[idx_s].set_offsets(x)
    # scatters[idx_s].set_UVC(plotables.transport_force[:,1], plotables.transport_force[:,2])
    # idx_s += 1

    if F.centrosome
        # Centrosome
        ## Orientation
        scatters[idx_s].set_offsets(coords.centro_x) # centrosome
        scatters[idx_s].set_UVC([cos.(coords.centro_angle)], [sin.(coords.centro_angle)])
        idx_s += 1

        ## Visibility
        patches[idx_p].set_xy(vr.nodes[1:vr.n,:] .+ reshape(coords.centro_x, 1, 2))
        idx_p += 1

        ## Center
        # artists[idx_a].center = (coords.centro_x[1], coords.centro_x[2])
        # idx_a += 1

        ## Microtubule force
        # scatters[idx_s].set_offsets(coords.centro_x) # centrosome
        # scatters[idx_s].set_UVC([plotables.mt_force[1]], [plotables.mt_force[2]])
        # idx_s += 1

        # scatters[idx_s].set_offsets(vr.nodes[1:vr.n,:].+reshape(coords.centro_x, 1, 2)) # MT force
        # scatters[idx_s].set_UVC(1e-1plotables.mt_force_indiv[1:vr.n,1], 1e-1plotables.mt_force_indiv[1:vr.n,2])
        # idx_s += 1
    end

    # Added mass
    # scatters[idx_s].set_offsets(x)
    # scatters[idx_s].set_sizes(0e4*plotables.mass_source)
    # colors = Array{String}(undef, size(x, 1))
    # fill!(colors, "red")
    # x_max, x_max_idx = findmax(x[:,2])
    # colors[x_max_idx] = "yellow"
    # scatters[idx_s].set_color(colors)
    # idx_s += 1


    if F.follow_cam
        if F.nucleus && F.follow_nucleus
            follow_x = nucleus_coords.Y
        else
            follow_x = coords.x
        end
        x_min, x_max = minimum(follow_x[:,1]), maximum(follow_x[:,1])
        y_min, y_max = minimum(follow_x[:,2]), maximum(follow_x[:,2])
        x_mid = 0.5(x_max+x_min)
        x_span = (x_max-x_min)
        y_mid = 0.5(y_max+y_min)
        y_span = (y_max-y_min)

        ax.set_xlim((x_mid-0.7*x_span, x_mid+0.7*x_span))
        ax.set_ylim((y_mid-0.7*y_span, y_mid+0.7*y_span))
    end

    PyPlot.draw()
    sleep(0.001)
end

function init_animation(prefix::String, date_string::String)
    FFMpegWriter = animation.FFMpegWriter
    metadata = Dict((:title => string(prefix, "Run_", date_string), :artist => "GJ"))
    writer = FFMpegWriter(fps=15, metadata=metadata)
    return writer
end

function plot_metrics(m::Dict)
    n_data_points = length(m["inst_max_y"])
    PyPlot.figure()
    PyPlot.subplot(231)
    println(typeof(m["inst_max_y"]))
    PyPlot.plot(m["inst_max_y"], label="tip")
    if haskey(m, "inst_bc")
        PyPlot.plot(m["inst_bc"][:,2], label="barycenter")
    end
    if haskey(m, "inst_cs")
        PyPlot.plot(m["inst_cs"][:,2], label="centrosome")
    end
    PyPlot.legend()
    PyPlot.title("Position")
    PyPlot.subplot(234)
    PyPlot.plot(m["inst_velocity"], label="tip")
    PyPlot.plot(m["inst_bc_velocity"], label="barycenter")
    PyPlot.title("Velocity")
    PyPlot.legend()
    PyPlot.subplot(233)
    PyPlot.plot(m["inst_cortex_area"], label="Cortex")
    if haskey(m, "inst_nucleus_area")
        PyPlot.plot(m["inst_nucleus_area"], label="Nucleus")
    end
    PyPlot.title("Area")
    PyPlot.legend()
    PyPlot.subplot(236)
    PyPlot.plot(4π*m["inst_cortex_area"]./(m["inst_cortex_perimeter"].^2), label="Cortex")
    if haskey(m, "inst_nucleus_area")
        PyPlot.plot(4π*m["inst_nucleus_area"]./(m["inst_nucleus_perimeter"].^2), label="nucleus")
    end
    PyPlot.axhline(1.0)
    PyPlot.title("Roundness (4π Area/Perimeter²)")
    PyPlot.legend()
    PyPlot.subplot(232)
    PyPlot.title("Distances")
    if haskey(m, "inst_n2c_distance")
        PyPlot.plot(m["inst_n2c_distance"], label="Nucleus barycenter to centrosome")
    end
    if haskey(m, "inst_bc2cs_distance")
        PyPlot.plot(m["inst_bc2cs_distance"], label="Cortex barycenter to centrosome")
    end
    PyPlot.legend()
    if PyPlot.matplotlib.backends.backend != "agg"
        PyPlot.show()
    end
end

function close()
    PyPlot.close("all")
end

end
