module CellSim

import CellSimCommon
import Forces
import Wall
import Masks

import Centrosome

import Utils
import SurfacesCommon
import EvenParam

import PyPlot
import PyCall
PyCall.@pyimport matplotlib.animation as animation

import NLsolve

import YAML

macro eval_if_string(s)
    return esc(:(isa($s, String) ? eval(parse($s)) : $s))
end

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

        ax[:scatter](x[:,1], x[:,2], color="black", zorder=30)
        ax[:plot](x[:,1], x[:,2], color="black", lw=0.5, zorder=1)[1] # initial condition

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

function init_animation()
    FFMpegWriter = animation.FFMpegWriter
    date = string(Dates.now())
    metadata = Dict((:title => string("Cell_run_", date), :artist => "GJ"))
    writer = FFMpegWriter(fps=15, metadata=metadata)
    return writer
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

    if !F.landscape_plot
        lines[1][:set_data](x[:,1], x[:,2])
        patches[1][:set_xy](x)

        if F.centrosome
            scatters[2][:set_offsets](coords.centro_x) # centrosome
            scatters[2][:set_UVC]([cos.(coords.centro_angle)], [sin.(coords.centro_angle)])
            patches[2][:set_xy](vr.nodes[1:vr.n,:] .+ reshape(coords.centro_x, 1, 2)) # visibility region

            scatters[3][:set_offsets](coords.centro_x) # centrosome
            scatters[3][:set_UVC]([plotables.mt_force[1]], [plotables.mt_force[2]])

            scatters[4][:set_offsets](vr.nodes[1:vr.n,:].+reshape(coords.centro_x, 1, 2)) # MT force
            scatters[4][:set_UVC]([plotables.mt_force_indiv[1:vr.n,1]], plotables.mt_force_indiv[1:vr.n,2])

            artists[1][:center] = (coords.centro_x[1], coords.centro_x[2])
        end

        # drag force
        scatters[1][:set_offsets](x[:,:])

        lims = ax[:get_ylim]()

        x_min, x_max = minimum(x[:,2]), maximum(x[:,2])
        x_mid = 0.5(x_max+x_min)
        # ax[:set_ylim]((x_mid-15, x_mid+15))
    else
        lines[1][:set_data](x[:,2], x[:,1])
        patches[1][:set_xy]([x[:,2] x[:,1]]) # cell polygon

        if F.centrosome
            scatters[2][:set_offsets]([coords.centro_x[2]; coords.centro_x[1]]) # centrosome
            scatters[2][:set_UVC]([sin.(coords.centro_angle)], [cos.(coords.centro_angle)])
            patches[2][:set_xy]([vr.nodes[1:vr.n,2]+coords.centro_x[2] vr.nodes[1:vr.n,1]+coords.centro_x[1]]) # visibility region

            scatters[3][:set_offsets]([coords.centro_x[2]], [coords.centro_x[1]]) # MT force
            scatters[3][:set_UVC]([plotables.mt_force[2]], [plotables.mt_force[1]])

            scatters[4][:set_offsets]([vr.nodes[1:vr.n,2]+coords.centro_x[2] vr.nodes[1:vr.n,1]+coords.centro_x[1]]) # MT force
            scatters[4][:set_UVC]([plotables.mt_force_indiv[1:vr.n,2]], plotables.mt_force_indiv[1:vr.n,1])

            artists[1][:center] = (coords.centro_x[2], coords.centro_x[1])
        end

        # drag force
        scatters[1][:set_offsets]([x[:,2] x[:,1]])

        lims = ax[:get_xlim]()
        x_min, x_max = minimum(x[:,2]), maximum(x[:,2])
        x_mid = 0.5(x_max+x_min)
        # ax[:set_xlim]((x_mid-15, x_mid+15))
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

function compute_initial_x(P::CellSimCommon.Params, F::CellSimCommon.Flags; convex::Bool=true)
    t = linspace(0, 1, P.N+1)[1:P.N]

    if convex
        if !F.circular_wall
            return 0.5 * Float64[P.x0_a*cospi.(2t) P.x0_b*sinpi.(2t)]
        else
            return 0.5 * Float64[P.x0_a*cospi.(2t)+2P.polar_shift P.x0_b*sinpi.(2t)]
        end
    end

    cx, cy = 0.6, 0.9
    xr = 0.24
    xl = 1-xr
    yr, yl = 0.25, 0.75
    pow = 2.0
    px, py = 1600, 600

    circ_x = cospi.(2t)
    bump_up = cx*exp.(-px*(t-xr).^pow)-cx*exp.(-px*(t-(0.5-xr)).^pow)
    bump_down =  cx*exp.(-px*(t-xl).^pow)-cx*exp.(-px*(t-(1.5-xl)).^pow)

    x = P.x0_a*(circ_x+bump_up+bump_down)
    y =  P.x0_b*(sinpi.(2t)-cy*exp.(-py*(t-yr).^pow)+cy*exp.(-py*(t-yl).^pow))

    if !F.circular_wall
        return Float64[x y]
    else
        return Float64[x+P.polar_shift y]
    end
end

function main()
    # initialization

    config_filename = "config.yaml"
    if length(ARGS) > 0 && isfile(ARGS[1])
        config_filename = ARGS[1]
    else
        println("Please provide a configuration file!")
        println("Usage: run.jl <config>")
        exit(0)
    end

    yaml_config = YAML.load(open(config_filename))
    y_params = yaml_config["params"]

    # Parameters
    P = CellSimCommon.Params(
        @eval_if_string(y_params["N"]),
        y_params["M"],
        1/y_params["N"], # Δσ
        y_params["δt"],
        @eval_if_string(y_params["P"]),
        @eval_if_string(y_params["K"]),
        @eval_if_string(y_params["Ka"]),
        @eval_if_string(y_params["c"]),
        y_params["x0_a"],
        y_params["x0_b"],
        y_params["x0_shift"],
        y_params["f_α"],
        y_params["f_β"],
        @eval_if_string(y_params["f_ω0"]),
        y_params["f_σ"],
        y_params["f_nk"],
        y_params["f_width"],
        y_params["f_iwidth"],
        y_params["drag_gauss_power"],
        y_params["drag_gauss_width"],
        y_params["mass_gauss_power"],
        y_params["mass_gauss_width"],
        y_params["polar_shift"],
        y_params["k_MT"]
   )

    println("equilibrium radius: ", 1/(2*pi - P.P/P.K))
    println("critical pressure: ", 2*pi*P.K)

    y_flags = yaml_config["flags"]

    # Flags
    F = CellSimCommon.Flags(
          # model options
          y_flags["confine"],
          y_flags["adjust_drag"],
          y_flags["polymerize"],
          y_flags["continuous"],
          y_flags["circular_wall"],
          y_flags["cortex"],
          y_flags["centrosome"],
          y_flags["weighted_confinement"],

          # scheme/solver options
          y_flags["innerloop"],

          # visualization options
          y_flags["plot"],
          y_flags["pretty"],
          y_flags["landscape_plot"],
          y_flags["follow_cam"],
          y_flags["plot_drag"],

          # output options
          y_flags["dryrun"],
          y_flags["write_animation"]
    )

    println(P)
    println(F)

    CellSimCommon.init(P.N)

    plotables = CellSimCommon.new_plotables(P.N)

    if haskey(yaml_config, "load_state") && yaml_config["load_state"]["do_load"]
        x_init = readcsv(yaml_config["load_state"]["filename"])
        x_init .-= sum(x_init, 1)/size(x_init, 1)
        x_init[:,1] *= P.x0_a/abs(maximum(x_init[:,1]))
        x_init[:,2] *= P.x0_b/abs(maximum(x_init[:,2]))

        if F.circular_wall
            x_init[:,1] += P.polar_shift
        end
        x_init = EvenParam.reparam(x_init, true, P.N)
        println("[info] reversing initial condition")
        if !Utils.check_ccw_polygon(x_init)
            @views begin
                reverse!(x_init[:,1])
                reverse!(x_init[:,2])
            end
        end
    else
        x_init = EvenParam.reparam(compute_initial_x(P, F; convex=true))
    end


    x = copy(x_init)

    Forces.init_FD_matrices(P)
    coords, coords_s = Forces.new_PointCoords(x, P)

    if haskey(yaml_config, "load_state") && yaml_config["load_state"]["init_centro"]
        # if the initial centrosome location is given in the config, load it
        # the centrosome angle doesn't have any impact at this point
        coords.centro_x[:] = [yaml_config["load_state"]["centro_x"]; yaml_config["load_state"]["centro_y"]]
    else
        # otherwise, pick the centrosome location as the initial center of mass
        coords.centro_x[:] = sum(x_init, 1)/size(x_init,1)
    end

    fig = init_plot(coords, P, F)
    if F.write_animation
        writer = init_animation()
        writer[:setup](fig, string(writer[:metadata]["title"], ".mp4"), 100)
    end

    # centrosome buffers and coordinates
    (centro_bufs, centro_vr, centro_qw, centro_pc) = Centrosome.init(P)

    resi, resi_J = Forces.wrap_residuals(coords, coords_s, P, F, plotables)
    if F.innerloop
        resi_solver = NLsolve.DifferentiableSparseMultivariateFunction(resi, resi_J)
    else
        r_x = zeros(2P.N)
        Jr_x = spzeros(2P.N,2P.N)
        δx = zeros(2P.N)
    end

    k = 0
    prev_height = 0.0

    # outer loop
    while k < P.M
        k += 1
        println("iteration #", k)

        # inner loop
        if k > 1
            Forces.update_coords(coords, P, x)
        end

        if F.innerloop
            res = NLsolve.nlsolve(resi_solver, vec(x); method=:newton)
            x = reshape(res.zero, (P.N, 2))
        else
            if F.cortex
                # cortex evolution
                resi(vec(x), r_x)
                resi_J(vec(x), Jr_x)
                δx[:] = -(Jr_x\r_x)
                x[:] = x[:] + δx
            end

            if F.centrosome
                # centrosome evolution
                (centro_A, centro_b) = Centrosome.assemble_system(P, coords, centro_bufs, centro_vr, centro_qw, centro_pc, plotables)
                centro_delta = centro_A\-centro_b
                coords.centro_x[:] += P.δt * centro_delta[1:2]
                coords.centro_angle[:] += P.δt * centro_delta[3]

                fill!(plotables.mt_force, 0.0)
                plotables.mt_force[:] = centro_b[1:2]
            end
        end


        # plot
        plot_period = 10
        if F.plot & (k % plot_period == 0)
            update_plot(coords, k, P, F, false, plotables, centro_vr)
            if F.write_animation
                writer[:grab_frame]()
            end
            height = sum(x[:,2]/P.N)
            long_speed = (height - prev_height) / (plot_period*P.δt)
            prev_height = height
            # println("Long. speed: ", long_speed)
        end

        l2_norm = sqrt(sum(abs2, x))

        # if l2_norm > 1e4
        # println("Divergence detected, aborting")
        # break
        # end
    end

    if F.write_animation
        writer[:finish]()
    end
    println("Finished, type Enter to exit")
    read(STDIN, 1)
    # output
end

end # module
