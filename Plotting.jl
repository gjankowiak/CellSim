module Plotting

import PyPlot
import PyCall
PyCall.@pyimport matplotlib.animation as animation

import CellSimCommon
import Centrosome
import Forces
import Wall
import Utils

function init_plot(coords::Forces.PointCoords, P::CellSimCommon.Params, F::CellSimCommon.Flags)
    PyPlot.ion()

    x = coords.x

    fig = PyPlot.figure(figsize=(12.8, 10))
    # figManager = PyPlot.get_current_fig_manager()
    # figManager[:window][:showMaximized]()

    if !F.landscape_plot
        ax = PyPlot.axes(xlim = (-5, 5),ylim=(-15, 15))

        ax[:set_aspect]("equal", "datalim")

        line = ax[:plot](x[:,1], x[:,2], ".-", zorder=20)[1]
        polygon = ax[:fill](x[:,1], x[:,2], color="#f713e0", zorder=10)[1]

        ax[:scatter](x[:,1], x[:,2], color="black", zorder=30) # ?
        ax[:plot](x[:,1], x[:,2], color="black", lw=0.5, zorder=1)[1] # initial condition

        ax[:quiver](coords.x[1], coords.x[2],
                    zeros(size(x,1)), zeros(size(x,1)), zorder=100, units="xy", scale=30, width=0.001, color="blue") # transport force

        if F.centrosome
            ax[:quiver](coords.centro_x[1], coords.centro_x[2],
                        [cos.(coords.centro_angle)], [sin.(coords.centro_angle)], zorder=100, units="xy", scale=10, width=0.01, color="#12b600") # centrosome
            ax[:fill](x[:,1], x[:,2], color="#32d600", zorder=15)[1] # centrosome visibility
            ax[:quiver](coords.centro_x[1], coords.centro_x[2],
                        [0], [0], zorder=100, units="xy", scale=1e-2, width=0.01) # MT force

            ax[:quiver]([0.0], [0.0], [0], [0], zorder=100, units="xy", scale=5, width=0.001) # MT force

            mt_circle = PyPlot.matplotlib[:patches][:Circle]((coords.centro_x[1], coords.centro_x[2]), Centrosome.MT_RADIUS, zorder=110, alpha=0.5)
            ax[:add_artist](mt_circle)
        end

        if F.circular_wall
            y = collect(linspace(-1, 1, 1000))
            wall_1 = Wall.compute_walls(y, P, F)
            wall_2 = Wall.compute_walls(y, P, F; right=true)
            ax[:plot](wall_1[:,1], wall_1[:,2], wall_2[:,1], wall_2[:,2], color="black", lw=0.5)

            levelset_1 = Wall.compute_walls(y, P, F, 2e-4)
            levelset_2 = Wall.compute_walls(y, P, F, 2e-4; right=true)
            ax[:plot](levelset_1[:,1], levelset_1[:,2], levelset_2[:,1], levelset_2[:,2], color="red", lw=0.5)
        else
            y = collect(linspace(-40, 40, 1000))
            wall = Wall.compute_walls(y, P, F)
            ax[:plot](wall[:,1], wall[:,2], -wall[:,1], wall[:,2], color="black", lw=0.5)
            levelset = Wall.compute_walls(y, P, F, 2e-4)
            ax[:plot](levelset[:,1], levelset[:,2], -levelset[:,1], levelset[:,2], color="red", lw=0.5)
        end

        ax[:scatter](x[:,1], x[:,2], color="red")

        ax[:axvline](0)
    else
        ax = PyPlot.axes(xlim = (-15, 15),ylim=(-5, 5))

        ax[:set_aspect]("equal", "datalim")

        line = ax[:plot](x[:,2], x[:,1], ".-", zorder=20)[1]
        polygon = ax[:fill](x[:,2], x[:,1], color="#f713e0", zorder=10)[1]

        ax[:scatter](x[:,2], x[:,1], color="black", zorder=30)
        # ax[:scatter](x[:,1], x[:,2], color="black", zorder=2)
        ax[:plot](x[:,2], x[:,1], color="black", lw=0.5)[1] # initial condition

        if F.centrosome
            ax[:quiver](coords.centro_x[2], coords.centro_x[1],
                        [sin.(coords.centro_angle)], [cos.(coords.centro_angle)], zorder=100, units="xy", scale=10, width=0.01, color="#12b600") # centrosome
            ax[:fill](x[:,2], x[:,1], color="#32d600", zorder=15)[1] # centrosome visibility
            ax[:quiver](coords.centro_x[2], coords.centro_x[1],
                        [0], [0], zorder=100, units="xy", scale=1e-2, width=0.01) # MT force

            ax[:quiver]([0.0], [0.0], [0], [0], zorder=100, units="xy", scale=5, width=0.001) # MT force

            mt_circle = PyPlot.matplotlib[:patches][:Circle]((coords.centro_x[2], coords.centro_x[1]), Centrosome.MT_RADIUS, zorder=110, alpha=0.5)
            ax[:add_artist](mt_circle)
        end

        if F.circular_wall
            y = collect(linspace(-1, 1, 1000))
            wall_1 = Wall.compute_walls(y, P, F)
            wall_2 = Wall.compute_walls(y, P, F; right=true)
            ax[:plot](wall_1[:,2], wall_1[:,1], wall_2[:,2], wall_2[:,1], color="black", lw=0.5)

            levelset_1 = Wall.compute_walls(y, P, F, 2e-4)
            levelset_2 = Wall.compute_walls(y, P, F, 2e-4; right=true)
            ax[:plot](levelset_1[:,2], levelset_1[:,1], levelset_2[:,2], levelset_2[:,1], color="red", lw=0.5)
        else
            y = collect(linspace(-40, 40, 1000))
            wall = Wall.compute_walls(y, P, F)
            ax[:plot](wall[:,1], wall[:,2], -wall[:,1], wall[:,2], color="black", lw=0.5)
            levelset = Wall.compute_walls(y, P, F, 2e-4)
            ax[:plot](levelset[:,2], levelset[:,1], -levelset[:,2], levelset[:,1], color="red", lw=0.5)
        end

        ax[:axhline](0)
    end

    lims = ax[:get_xlim]()
    x_min, x_max = minimum(x[:,1]), maximum(x[:,1])
    x_mid = 0.5(x_max+x_min)
    x_span = (x_max-x_min)
    y_min, y_max = minimum(x[:,2]), maximum(x[:,2])
    y_mid = 0.5(y_max+y_min)
    y_span = (y_max+y_min)

    if !F.landscape_plot
        ax[:set_xlim]((x_mid-x_span, x_mid+x_span))
        ax[:set_ylim]((y_mid-y_span, y_mid+y_span))
    else
        ax[:set_ylim]((x_mid-x_span, x_mid+x_span))
        ax[:set_xlim]((y_mid-y_span, y_mid+y_span))
    end

    #PyPlot.show()

    fig[:tight_layout]()
    println("Finishing plot init...")
    sleep(1)
    return fig
end

function update_plot(coords::Forces.PointCoords, k::Int, P::CellSimCommon.Params, F::CellSimCommon.Flags, initializing::Bool, plotables::CellSimCommon.Plotables,
                     vr::Centrosome.VisibleRegion)

    x = coords.x

    ax = PyPlot.gca()
    if initializing
        prefix = "[INIT]"
    else
        prefix = ""
    end
    ax[:set_title](@sprintf("%s N: %d, iter: %d, T=%fs", prefix, P.N, k, k*P.δt))

    lines = ax[:lines]
    scatters = ax[:collections]
    patches = ax[:patches]
    artists = ax[:artists]

    idx_l = idx_s = idx_p = idx_a = 1

    if !F.landscape_plot
        lines[idx_l][:set_data](x[:,1], x[:,2])
        idx_l += 1

        patches[idx_p][:set_xy](x)
        idx_p += 1

        scatters[idx_s][:set_offsets](x)
        idx_s += 1

        scatters[idx_s][:set_offsets](x)
        scatters[idx_s][:set_UVC](plotables.transport_force[:,1], plotables.transport_force[:,2])
        idx_s += 1

        if F.centrosome
            scatters[idx_s][:set_offsets](coords.centro_x) # centrosome
            scatters[idx_s][:set_UVC]([cos.(coords.centro_angle)], [sin.(coords.centro_angle)])
            idx_s += 1

            patches[idx_p][:set_xy](vr.nodes[1:vr.n,:] .+ reshape(coords.centro_x, 1, 2)) # visibility region
            idx_p += 1

            scatters[idx_s][:set_offsets](coords.centro_x) # centrosome
            scatters[idx_s][:set_UVC]([plotables.mt_force[1]], [plotables.mt_force[2]])
            idx_s += 1

            scatters[idx_s][:set_offsets](vr.nodes[1:vr.n,:].+reshape(coords.centro_x, 1, 2)) # MT force
            scatters[idx_s][:set_UVC]([plotables.mt_force_indiv[1:vr.n,1]], plotables.mt_force_indiv[1:vr.n,2])
            idx_s += 1

            artists[idx_a][:center] = (coords.centro_x[1], coords.centro_x[2])
            idx_a += 1
        end

        scatters[idx_s][:set_offsets](x)
        scatters[idx_s][:set_sizes](1e4*plotables.mass_source)
        colors = Array{String}(size(x, 1))
        fill!(colors, "red")
        x_max, x_max_idx = findmax(x[:,2])
        colors[x_max_idx] = "yellow"
        scatters[idx_s][:set_color](colors)
        idx_s += 1

        lims = ax[:get_ylim]()

        x_min, x_max = minimum(x[:,2]), maximum(x[:,2])
        x_mid = 0.5(x_max+x_min)
    else
        lines[idx_l][:set_data](x[:,2], x[:,1])
        idx_l += 1

        patches[idx_p][:set_xy]([x[:,2] x[:,1]]) # cell polygon
        idx_p += 1

        if F.centrosome
            scatters[idx_s][:set_offsets]([coords.centro_x[2]; coords.centro_x[1]]) # centrosome
            scatters[idx_s][:set_UVC]([sin.(coords.centro_angle)], [cos.(coords.centro_angle)])
            idx_s += 1

            patches[idx_p][:set_xy]([vr.nodes[1:vr.n,2]+coords.centro_x[2] vr.nodes[1:vr.n,1]+coords.centro_x[1]]) # visibility region
            idx_p += 1

            scatters[idx_s][:set_offsets]([coords.centro_x[2]], [coords.centro_x[1]]) # MT force
            scatters[idx_s][:set_UVC]([plotables.mt_force[2]], [plotables.mt_force[1]])
            idx_s += 1

            scatters[idx_s][:set_offsets]([vr.nodes[1:vr.n,2]+coords.centro_x[2] vr.nodes[1:vr.n,1]+coords.centro_x[1]]) # MT force
            scatters[idx_s][:set_UVC]([plotables.mt_force_indiv[1:vr.n,2]], plotables.mt_force_indiv[1:vr.n,1])
            idx_s += 1

            artists[idx_a][:center] = (coords.centro_x[2], coords.centro_x[1])
            idx_a += 1
        end

        lims = ax[:get_xlim]()
        x_min, x_max = minimum(x[:,2]), maximum(x[:,2])
        x_mid = 0.5(x_max+x_min)
    end


    if F.follow_cam
        lims = ax[:get_xlim]()
        x_min, x_max = minimum(x[:,1]), maximum(x[:,1])
        x_mid = 0.5(x_max+x_min)
        x_span = (x_max-x_min)
        y_min, y_max = minimum(x[:,2]), maximum(x[:,2])
        y_mid = 0.5(y_max+y_min)
        y_span = (y_max-y_min)

        if !F.landscape_plot
            ax[:set_xlim]((x_mid-x_span, x_mid+x_span))
            ax[:set_ylim]((y_mid-y_span, y_mid+y_span))
        else
            ax[:set_ylim]((x_mid-x_span, x_mid+x_span))
            ax[:set_xlim]((y_mid-y_span, y_mid+y_span))
        end
    end

    if F.plot_drag
        scatters[1][:set_sizes](10CellSimCommon.@entry_norm(plotables.drag_force))
    else
        scatters[1][:set_sizes](0CellSimCommon.@entry_norm(plotables.drag_force))
    end
    scatters[1][:set_facecolor](Utils.scale_cm(plotables.mass_source, PyPlot.get_cmap("RdYlGn");
                                               range_min=Nullable(-P.c), range_max=Nullable(P.c)))


    PyPlot.draw()
    sleep(0.0001)
end

function init_animation()
    FFMpegWriter = animation.FFMpegWriter
    date = string(Dates.now())
    metadata = Dict((:title => string("Cell_run_", date), :artist => "GJ"))
    writer = FFMpegWriter(fps=15, metadata=metadata)
    return writer
end


end
