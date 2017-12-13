module Refactor

import AdhCommon
import Forces
import Wall
import Masks
import Utils
# import Surface
import SurfacesCommon
import EvenParam

import PyPlot
import NLsolve

import YAML

macro eval_if_string(s)
    return esc(:(isa($s, String) ? eval(parse($s)) : $s))
end

function init_plot(x::Matrix, P::AdhCommon.Params, F::AdhCommon.Flags)
    PyPlot.ion()

    fig = PyPlot.figure()
    ax = PyPlot.axes(xlim = (-10, 10),ylim=(-15, 15))

    ax[:set_aspect]("equal", "datalim")

    line = ax[:plot](x[:,1], x[:,2], ".-", zorder=1)[1]
    ax[:scatter](x[:,1], x[:,2], color="black", zorder=2)
    # ax[:scatter](x[:,1], x[:,2], color="black", zorder=2)
    ax[:plot](x[:,1], x[:,2], color="black", lw=0.5)[1] # initial condition

    y = collect(linspace(-40, 40, 1000))
    wall = Wall.compute_walls(y, P)
    levelset = Wall.compute_walls(y, P, 1e-8)
    ax[:plot](wall, y, -wall, y, color="black", lw=0.5)

    levelset = Wall.compute_walls(y, P, 2e-4)
    ax[:plot](levelset, y, -levelset, y, color="red", lw=0.5)

    ax[:axvline](0)

    # figManager = PyPlot.get_current_fig_manager()
    # figManager[:window][:showMaximized]()
    PyPlot.show()

    fig[:tight_layout]()
    return line
end

function update_plot(x::Matrix, k::Int, P::AdhCommon.Params, F::AdhCommon.Flags, initializing::Bool, plotables::Forces.Plotables)
    ax = PyPlot.gca()
    lines = ax[:lines]
    scatters = ax[:collections]
    lines[1][:set_data](x[:,1], x[:,2])
    # tracker[:set_data](x[tr_idx,1], x[tr_idx,2])

    # line_dbg[:set_data](1:N, x[:,1]-x0[:,1])

    # effective angle force
    # FA = @entry_norm(reshape(M_DFangle * δx, N, 2))/δt

    if initializing
        prefix = "[INIT]"
    else
        prefix = ""
    end
    ax[:set_title](@sprintf("%s N: %d, iter: %d, T=%fs", prefix, P.N, k, k*P.δt))

    # effective angle force
    # FA = @entry_norm(reshape(M_DFangle * δx, N, 2))/δt

    # drag force
    scatters[1][:set_offsets](x[:,:])
    scatters[1][:set_sizes](20AdhCommon.@entry_norm(plotables.drag_force))
    scatters[1][:set_facecolor](Utils.scale_cm(plotables.mass_source, PyPlot.get_cmap("RdYlGn");
                                                range_min=Nullable(-P.c), range_max=Nullable(P.c)))

    # mass source
    # center = sum(x, 1)/P.N
    # offset_from_center = broadcast(-, x, center)
    # scatters[2][:set_offsets](broadcast(+, center, 0.95offset_from_center))
    # scatters[2][:set_sizes](50abs.(plotables.mass_source))

    PyPlot.draw()
    sleep(0.0001)
    # read(STDIN, 1)
end

function main()
    # initialization
    yaml_config = YAML.load(open("config.yaml"))
    y_params = yaml_config["params"]

    # Parameters
    P = AdhCommon.Params(
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
            y_params["f_ω0"],
            y_params["f_σ"],
            y_params["f_nk"],
            y_params["f_width"],
            y_params["f_iwidth"],
            y_params["drag_gauss_power"],
            y_params["drag_gauss_width"],
            y_params["mass_gauss_power"],
            y_params["mass_gauss_width"]
            )

    println("equilibrium radius: ", 1/(2*pi - P.P/P.K))
    println("critical pressure: ", 2*pi*P.K)

    y_flags = yaml_config["flags"]

    # Flags
    F = AdhCommon.Flags(
            y_flags["confine"],
            y_flags["adjust_drag"],
            y_flags["polymerize"],
            y_flags["dryrun"],
            y_flags["plot"],
            y_flags["pretty"],
            y_flags["continuous"],
            y_flags["innerloop"],
            y_flags["weighted_confinement"]
           )

    println(P)
    println(F)

    AdhCommon.init(P.N)

    plotables = Forces.new_plotables(P.N)

    t = linspace(0, 1, P.N+1)[1:P.N]

    x_init = 0.5 * Float64[P.x0_a*cospi.(2t) P.x0_b*sinpi.(2t)]
    x = EvenParam.reparam(x_init, true)
    line = init_plot(x_init, P, F)

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
        plot_period = 20
        if F.plot & (k % plot_period == 0)
            update_plot(x, k, P, F, false, plotables)
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

    println("Finished, type Enter to exit")
    read(STDIN, 1)
    # output
end

end # module
