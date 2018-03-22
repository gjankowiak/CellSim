module CellSim

import CellSimCommon
import Forces
import Wall
import Masks
import Utils
# import Surface
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

function init_plot(x::Matrix, P::CellSimCommon.Params, F::CellSimCommon.Flags)
    PyPlot.ion()

    fig = PyPlot.figure(figsize=(12.8, 10))
    # figManager = PyPlot.get_current_fig_manager()
    # figManager[:window][:showMaximized]()

    if !F.landscape_plot
        ax = PyPlot.axes(xlim = (-5, 5),ylim=(-15, 15))

        ax[:set_aspect]("equal", "datalim")

        line = ax[:plot](x[:,1], x[:,2], ".-", zorder=1)[1]
        polygon = ax[:fill](x[:,2], x[:,1], color="#f713e0", zorder=1)[1]
        ax[:scatter](x[:,1], x[:,2], color="black", zorder=2)
        # ax[:scatter](x[:,1], x[:,2], color="black", zorder=2)
        ax[:plot](x[:,1], x[:,2], color="black", lw=0.5)[1] # initial condition



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
            ax[:plot](levelset[:,1], levelset[:,2], -levelset[:,1], levelset[:,1], color="red", lw=0.5)
        end

        ax[:axvline](0)
    else
        ax = PyPlot.axes(xlim = (-15, 15),ylim=(-5, 5))

        ax[:set_aspect]("equal", "datalim")

        line = ax[:plot](x[:,2], x[:,1], ".-", zorder=1)[1]
        polygon = ax[:fill](x[:,2], x[:,1], color="#f713e0", zorder=1)[1]
        ax[:scatter](x[:,2], x[:,1], color="black", zorder=2)
        # ax[:scatter](x[:,1], x[:,2], color="black", zorder=2)
        ax[:plot](x[:,2], x[:,1], color="black", lw=0.5)[1] # initial condition

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

function update_plot(x::Matrix, k::Int, P::CellSimCommon.Params, F::CellSimCommon.Flags, initializing::Bool, plotables::Forces.Plotables)
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

    if !F.landscape_plot
        lines[1][:set_data](x[:,1], x[:,2])
        patches[1][:set_xy](x)

        # drag force
        scatters[1][:set_offsets](x[:,:])

        lims = ax[:get_ylim]()

        x_min, x_max = minimum(x[:,2]), maximum(x[:,2])
        x_mid = 0.5(x_max+x_min)
        # ax[:set_ylim]((x_mid-15, x_mid+15))
    else
        lines[1][:set_data](x[:,2], x[:,1])
        patches[1][:set_xy]([x[:,2] x[:,1]])

        # drag force
        scatters[1][:set_offsets]([x[:,2] x[:,1]])

        lims = ax[:get_xlim]()
        x_min, x_max = minimum(x[:,2]), maximum(x[:,2])
        x_mid = 0.5(x_max+x_min)
        # ax[:set_xlim]((x_mid-15, x_mid+15))
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
          y_flags["confine"],
          y_flags["adjust_drag"],
          y_flags["polymerize"],
          y_flags["dryrun"],
          y_flags["plot"],
          y_flags["pretty"],
          y_flags["continuous"],
          y_flags["innerloop"],
          y_flags["weighted_confinement"],
          y_flags["write_animation"],
          y_flags["landscape_plot"],
          y_flags["plot_drag"],
          y_flags["circular_wall"],
          y_flags["microtubules"]
    )

    println(P)
    println(F)

    CellSimCommon.init(P.N)

    plotables = Forces.new_plotables(P.N)

    x_init = EvenParam.reparam(compute_initial_x(P, F; convex=true))

    x = copy(x_init)
    fig = init_plot(x_init, P, F)
    if F.write_animation
        writer = init_animation()
        writer[:setup](fig, string(writer[:metadata]["title"], ".mp4"), 100)
    end

    Forces.init_FD_matrices(P)
    coords, coords_s = Forces.new_PointCoords(x, P)

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
        # println("iteration #", k)

        # inner loop
        if k > 1
            Forces.update_coords(coords, P, x)
        end

        if F.innerloop
            res = NLsolve.nlsolve(resi_solver, vec(x); method=:newton)
            x = reshape(res.zero, (P.N, 2))
        else
            resi(vec(x), r_x)
            resi_J(vec(x), Jr_x)
            δx[:] = -(Jr_x\r_x)
            x[:] = x[:] + δx
        end


        # plot
        plot_period = 10
        if F.plot & (k % plot_period == 0)
            update_plot(x, k, P, F, false, plotables)
            if F.write_animation
                writer[:grab_frame]()
            end
            height = sum(x[:,2]/P.N)
            long_speed = (height - prev_height) / (plot_period*P.δt)
            prev_height = height
            println("Long. speed: ", long_speed)
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
