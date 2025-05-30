module CellSim

import CellSimCommon
const CSC = CellSimCommon

import Cortex
import Nucleus
import Wall
import Masks

import Metrics

import Plotting

import Centrosome

import LinearAlgebra
const LA = LinearAlgebra

import Distributed

import JankoUtils
import EvenParam

import NLsolve

import YAML
import Dates

# we need python's YAML to dump config files
import PyCall
const pyyaml = PyCall.pyimport("yaml")

import DelimitedFiles: readdlm, writedlm

import SparseArrays
const SA = SparseArrays

macro eval_if_string(s)
    return esc(:(isa($s, String) ? eval(Meta.parse($s)) : $s))
end

function compute_initial_x(P::CSC.Params, F::CSC.Flags; fill_wall::Bool=false)
    t = collect(range(0; stop=1, length=P.N+1))[1:P.N]

    if !fill_wall
        if !F.circular_wall
            return (0.5 * Float64[P.x0_a*cospi.(2t) P.x0_b*sinpi.(2t)], [0.0; 0.0])
        else
            return (0.5 * Float64[P.x0_a*cospi.(2t) .+ 2P.polar_shift P.x0_b*sinpi.(2t)], [0.0; 0.0])
        end
    else
        # this should probably move to Walls.jl
        init_width = P.f_width - 0.9/P.f_α
        @assert init_width > 0
        f = (x::Array{Float64,1}) -> [P.target_area - 2*((init_width)*x[1] + P.f_β/P.f_ω0*sin(P.f_ω0*x[1] - 0.5pi))]
        fp = (x::Array{Float64,1}) -> [-2*((init_width) - P.f_β*cos(P.f_ω0*x[1]-0.5pi))]
        res = NLsolve.nlsolve(f, fp, [P.target_area/(2*(init_width))]; method=:broyden)
        init_min_y = -0.5pi/P.f_ω0
        init_max_y = res.zero[1] + init_min_y
        width_at_front = (init_width) + P.f_β*sin(P.f_ω0*init_max_y)
        @assert width_at_front > 0
        approx_cell_perimeter = init_width - P.f_β + width_at_front + 2*init_max_y
        np_side = Int(round(P.N * init_max_y/approx_cell_perimeter))
        np_front = Int(round(P.N * width_at_front/approx_cell_perimeter))
        np_back = P.N - 2*np_side - np_front
        if np_back < 4
            np_front = div(P.N - 2*np_side, 2)
            np_back = P.N - 2*np_side - np_front
        end

        y = collect(range(init_min_y, stop=init_max_y, length=np_side))
        x_init = [[init_width.+P.f_β*sin.(P.f_ω0*y) y];
                  [range(width_at_front, stop=-width_at_front, length=np_front) init_max_y*ones(np_front)];
                  [-init_width.-P.f_β*sin.(P.f_ω0*reverse(y)) reverse(y)];
                  [range(-init_width + P.f_β, stop=init_width - P.f_β, length=np_back) init_min_y*ones(np_back)]]
        mid_pulse = 2π/P.f_ω0*(0.25+round(0.5*res.zero[1]/(2π/P.f_ω0) - 0.5, RoundNearest))
        return (x_init, [0.0; mid_pulse])
    end
end

function main()
    config_filename = "configs/default.yaml"
    if length(ARGS) > 0 && isfile(ARGS[1])
        config_filename = ARGS[1]
    else
        println("Please provide a configuration file!")
        println("Usage: run.jl <config>")
        exit(0)
    end

    P, F, config = read_config(config_filename)
    launch(P, F, config)
end

function read_config(config_filename::String)
    yaml_config = YAML.load(open(config_filename))
    y_params = yaml_config["params"]

    # Parameters
    P = CSC.Params(
        @eval_if_string(y_params["M"]),
        1/@eval_if_string(y_params["N"]), # Δσ
        @eval_if_string(y_params["δt"]),

        # Cortex related parameters
        @eval_if_string(y_params["N"]),
        @eval_if_string(y_params["P"]),
        @eval_if_string(y_params["K"]),
        @eval_if_string(y_params["Ka"]),
        @eval_if_string(y_params["c"]),
        @eval_if_string(y_params["mu"]),
        @eval_if_string(y_params["target_area"]),

        # Initial condition parameters
        @eval_if_string(y_params["x0_a"]),
        @eval_if_string(y_params["x0_b"]),
        @eval_if_string(y_params["x0_shift"]),

        # Confinement field
        @eval_if_string(y_params["f_α"]),
        @eval_if_string(y_params["f_β"]),
        @eval_if_string(y_params["f_ω0"]),
        @eval_if_string(y_params["f_σ"]),
        @eval_if_string(y_params["f_nk"]),
        @eval_if_string(y_params["f_width"]),
        @eval_if_string(y_params["f_iwidth"]),

        @eval_if_string(y_params["drag_gauss_power"]),
        @eval_if_string(y_params["drag_gauss_width"]),

        @eval_if_string(y_params["mass_gauss_power"]),
        @eval_if_string(y_params["mass_gauss_width"]),

        @eval_if_string(y_params["polar_shift"]),

        # Centrosome related parameters
        @eval_if_string(y_params["k_MT"]),
        @eval_if_string(y_params["MT_potential_power"]),
        @eval_if_string(y_params["MT_factor"]),

        # Nucleus related parameters
        @eval_if_string(y_params["Nnuc"]),
        @eval_if_string(y_params["N_P"]),
        @eval_if_string(y_params["N_mu"]),
        @eval_if_string(y_params["N_target_area"]),
        @eval_if_string(y_params["N_kb"]),
        @eval_if_string(y_params["N_ω"]),
        @eval_if_string(y_params["N_W0"]),
        @eval_if_string(y_params["N_kcont"]),
        @eval_if_string(y_params["N_αcont"]),
        @eval_if_string(y_params["N_kc"]),
        @eval_if_string(y_params["N_l0c"]),
        @eval_if_string(y_params["N_r_init"])
   )

    y_flags = yaml_config["flags"]

    # Flags
    F = CSC.Flags(
          # model options
          y_flags["confine"],
          y_flags["adjust_drag"],
          y_flags["polymerize"],
          y_flags["continuous"],
          y_flags["circular_wall"],
          y_flags["cortex"],
          y_flags["centrosome"],
          y_flags["nucleus"],
          y_flags["weighted_confinement"],
          y_flags["force_cortex_area"],
          get(y_flags, "explicit_mt_force", false),

          # scheme/solver options
          y_flags["innerloop"],

          # visualization options
          y_flags["plot"],
          y_flags["pretty"],
          y_flags["landscape_plot"],
          y_flags["follow_cam"],
          get(y_flags, "follow_nucleus", false),
          y_flags["plot_drag"],

          # output options
          y_flags["dryrun"],
          y_flags["write_animation"],
          y_flags["write_metrics"],

          # debug
          y_flags["debug"]
    )

    config = yaml_config["config"]
    config["config_filename"] = config_filename

    if haskey(config, "output_prefix")
        if !endswith(config["output_prefix"], "/")
            config["output_prefix"] = config["output_prefix"] * "/"
        end
    else
        config["output_prefix"] = "runs/"
    end

    if !haskey(config, "recompute_nucleus_each")
        config["recompute_nucleus_each"] = 5
    end

    if (config["metrics"]["start_iteration"] <= 0) && !F.nucleus
        println("[WARNING] metrics.start_iteration is non-positive and nucleus is deactivated, metrics will never start!")
    end

    return P, F, config
end

function launch(P::CSC.Params, F::CSC.Flags, config; force_date_string::String="")

    if length("") > 0
        date_string = force_date_string
    else
        date_string = string(Dates.now())
    end

    config["date_string"] = date_string

    if haskey(config, "output_prefix")
        if !endswith(config["output_prefix"], "/")
            config["output_prefix"] = config["output_prefix"] * "/"
        end
    else
        config["output_prefix"] = "runs/"
    end

    mkpath("runs")
    mkpath("dumps")
    mkpath("$(config["output_prefix"])")
    config_dict = Dict([("params", CSC.to_dict(P)), ("flags", CSC.to_dict(F)), ("config", config)])
    output_config_file = open(config["output_prefix"] * "Run_" * config["date_string"] * ".yaml", "w")
    pyyaml.dump(config_dict, output_config_file, allow_unicode=true)
    close(output_config_file)
    println("Parameters:")
    for s in fieldnames(typeof(P))
        println(string("    ", s, ": ", getfield(P, s)))
    end
    println("")

    println("Flags:")
    for s in fieldnames(typeof(F))
        println(string("    ", s, ": ", getfield(F, s)))
    end
    println(" ")

    if F.debug
        println("*** debug MODE ***")
    end

    distr = (Distributed.nprocs() > 1)

    if !distr
        run(`stty -icanon`)
    end

    plotables = CSC.new_plotables(P.N)

    if haskey(config, "load_state")
        load_state = config["load_state"]
        if config["load_state"]["do_load"]
            x_init = readdlm(load_state["filename"], ',')

            if load_state["do_resample"]
                x_init = EvenParam.reparam(x_init; closed=true, new_N=P.N)
            end

            if !JankoUtils.check_ccw_polygon(x_init)
                println("[info] reversing initial condition")
                @views begin
                    reverse!(x_init[:,1])
                    reverse!(x_init[:,2])
                end
            end
        else

            (x_init, centro_x_init) = compute_initial_x(P, F; fill_wall=true)
            x_init[:] = EvenParam.reparam(x_init)
            x_init[:,2] .+= P.x0_shift
        end
    end


    x = copy(x_init)
    x_candidate = copy(x_init)

    solver_params = CSC.SolverParams(
            P.δt, # current time step
            1e-3, # maximum time step
            1e-6, # minimum time step
            1.4,  # stepping factor
            1e-4,  # step up error
            8e-4  # step down error
    )

    Cortex.init_FD_matrices(P)
    coords, coords_s = Cortex.new_PointCoords(x, P)
    old_coords = deepcopy(coords)

    if haskey(config, "load_state") && load_state["do_load"] && load_state["init_centro"]
        # if the initial centrosome location is given in the config, load it
        # the centrosome angle doesn't have any impact at this point
        coords.centro_x[:] = readdlm(load_state["filename_centro"], ',')
    else
        if F.nucleus
            # initizalize in the center of the nucleus
            coords.centro_x[:] = [0.0 0.5π/P.f_ω0]
        else
            # otherwise, pick the centrosome location as the initial center of mass
            coords.centro_x[:] = sum(x_init; dims=1)/size(x_init,1)
            if !F.circular_wall
                coords.centro_x[:] = centro_x_init
            end
        end
    end

    nucleus_coords = missing

    if F.nucleus
        nucleus_coords = Nucleus.initialize_coords(P, F, coords; fill_wall=true)
        if haskey(config, "load_state") && haskey(load_state, "shift_nucleus")
            nucleus_coords.Y[:] = nucleus_coords.Y .- sum(nucleus_coords.Y; dims=1)/P.Nnuc .+ [0.0 load_state["shift_nucleus"]]
        end
        old_nucleus_coords = deepcopy(nucleus_coords)
        #
        # shift nucleus to a wide section
        #

        temparrays = CSC.TempArrays6(
                                     zeros(P.Nnuc),
                                     zeros(P.Nnuc),
                                     zeros(P.Nnuc),
                                     zeros(P.Nnuc),
                                     zeros(P.Nnuc),
                                     zeros(P.Nnuc)
                                    )
    end

    potentials = CSC.InteractionPotentials(
        zeros(P.Nnuc),
        zeros(P.Nnuc, 2),
        zeros(P.N, 2),
        zeros(1, 2)
   )

    if F.centrosome
        # centrosome buffers and coordinates
        (centro_bufs, centro_vr, centro_qw, centro_pc) = Centrosome.init(P)
        Centrosome.compute_vr(P, coords, centro_bufs, centro_vr)
        (centro_A, centro_id_comp, centro_b_ce, centro_b_ce_rhs, centro_b_co_rhs) = Centrosome.assemble_system(P, F, coords, centro_bufs, centro_vr, centro_qw, centro_pc, plotables, potentials)
    else
        centro_vr = missing
    end

    centro_x_candidate = copy(coords.centro_x)
    centro_angle_candidate = copy(coords.centro_angle)

    resi, resi_J = Cortex.wrap_residuals(coords, coords_s, potentials, P, F, plotables)
    if F.innerloop
        resi_solver = NLsolve.DifferentiableSparseMultivariateFunction(resi, resi_J)
    else
        r_x = zeros(2P.N)
        Jr_x = SA.spzeros(2P.N,2P.N)
        if !F.centrosome
            δx = zeros(2P.N)
        else
            δx = zeros(2P.N+3)
        end
    end

    if F.plot
        if haskey(config, "plot_period")
            plot_period = config["plot_period"]
        end
        plot_period = F.debug ? 1 : 10
        fig = Plotting.init_plot(coords, P, F)
        Plotting.update_plot(coords, nucleus_coords, 0, 0.0, P, F, false, plotables, centro_vr)
        if F.write_animation
            writer = Plotting.init_animation(config["output_prefix"], config["date_string"])
            writer.setup(fig, string(writer.metadata["title"], ".mp4"), 100)
        end
    end

    k = 0
    prev_height = 0.0
    current_time = 0

    metrics = Metrics.init_metrics(P, F, config)
    post_init_periods = get(config["metrics"], "post_init_periods", 0)
    metrics_pre_init_done = false
    post_init_target_max_y = Inf

    stepping = F.debug
    breakpoint = -1
    if !stepping && !distr
        input_task = @async read(stdin, Char)
    end

    println("Initialized")

    # outer loop
    while k < P.M
        max_y = maximum(coords.x[:,2])
        if metrics["started"] && (max_y >= metrics["target_max_y"])
            metrics["finished"] = true
            println()
            print("Target number of periods reached, finishing...")
            println()
            break
        end

        if k == breakpoint
            stepping = true
        end

        if !distr
            if stepping
                println("debug: 'q' to quit, 'c' to run continuously, 'C' to continue to ITER, any other key to step, 'b' to run step by step")
                key = read(stdin, 1)[1]
                if key == 0x63 # c
                    stepping = false
                    input_task = @async read(stdin, Char)
                elseif key == 0x43 # C
                    stepping = false
                    breakpoint = parse(Int, readline())
                    input_task = @async read(stdin, Char)
                elseif key == 0x71 # q
                    break
                end
            elseif istaskdone(input_task)
                key = fetch(input_task)
                if key == 'b'
                    stepping = true
                elseif key == 'q'
                    break
                end
                input_task = @async read(stdin, Char)
            end
        end

        k += 1
        print("\b"^100)
        print(" iteration #", k, ", ")

        # Metrics initialization
        #
        # There are 2 ways to initialize:
        #    - using a fixed starting iteration. This uses the key "start_iteration" in the config file.
        #      If the value is positive, it will be the method used.
        #    - using a 2 stage initialization based on the distance between cortex and nucleus
        #      We detect the time when the distance between the back of the nucleus and the back
        #      of the cortex is smaller that the interaction distance (2/P.f_α)
        #
        #   Once one of previous condition is met, we then wait for the cell to cross
        #   "post_init_periods" after which we start recording metrics for "periods" periods.
        #
        if !metrics["started"]
                if !metrics_pre_init_done
                    if config["metrics"]["start_iteration"] > 0
                        if k == config["metrics"]["start_iteration"]
                            metrics_pre_init_done = true
                        end
                    else
                        if (F.nucleus && abs(minimum(coords.x[:,2]) - minimum(nucleus_coords.Y[:,2])) < 2/P.f_α)
                            metrics_pre_init_done = true
                        end
                    end
                    if metrics_pre_init_done
                        post_init_target_max_y = maximum(coords.x[:,2]) + 2π/P.f_ω0
                        println()
                        print("Starting post initialization, starting when head reaches y=", post_init_target_max_y)
                        println()
                    end
                elseif maximum(coords.x[:,2]) >= post_init_target_max_y
                    Metrics.start_metrics!(metrics, k, current_time, P, F, config, coords, nucleus_coords)
                    println()
                    print("Starting collecting metrics, stopping when head reaches y=", metrics["target_max_y"], " (+", config["metrics"]["periods"], ")")
                    println()
                elseif current_time > 3 # wait for max 20k iteration if P.M = 80k
                    println()
                    println("Initialization not finished after " * string(current_time) * " time units, aborting this run")
                    break
                end
        end


        try
        # inner loop
        if k > 1
            Cortex.update_coords(coords, P, x)
        end

        decreasing_step = false

        while true
            # adaptive

            if F.nucleus
                fill!(potentials.N_W, 0.0)
                fill!(potentials.N_∇W, 0.0)
                fill!(potentials.C_∇W, 0.0)
                fill!(potentials.CS_∇W, 0.0)
                Nucleus.compute_contact_force(potentials, coords, nucleus_coords, P, F)
                if F.centrosome
                    Nucleus.compute_centronuclear_force(potentials, coords, nucleus_coords, P, F)
                end
                recompute = (k % config["recompute_nucleus_each"] == 0)
                Nucleus.update_coords(old_nucleus_coords, nucleus_coords, potentials, P, F, temparrays, recompute)
            end

            if F.cortex
                # cortex evolution
                resi(vec(x), r_x)
                resi_J(vec(x), Jr_x)
            end

            if F.centrosome
                (centro_A, centro_id_comp, centro_b_ce, centro_b_ce_rhs, centro_b_co_rhs) = Centrosome.assemble_system(P, F, coords, centro_bufs, centro_vr, centro_qw, centro_pc, plotables, potentials)
                # centrosome evolution
                M = ([[Jr_x+centro_id_comp centro_b_ce'];[centro_b_ce centro_A]])

                rhs = [-r_x+P.δt*centro_b_co_rhs; P.δt*centro_b_ce_rhs]

                δx[:] = M\rhs

                x_candidate .= x .+ reshape(δx[1:2P.N], P.N, 2)

                centro_x_candidate .= coords.centro_x .+ δx[(2P.N+1):(2P.N+2)]
                centro_angle_candidate .= coords.centro_angle .+ δx[2P.N+3]
            else
                δx[:] = -Jr_x\r_x
                x_candidate .= x .+ reshape(δx, P.N, 2)
            end

            max_displacement = maximum(abs.(x - x_candidate))

            accept = true

            if max_displacement < solver_params.step_up_error # too small
                if decreasing_step
                    accept = true
                else
                    P.δt = min(solver_params.max_δt, P.δt * solver_params.stepping_factor)
                    # println()
                    # println("Time step ↑ $(P.δt)")
                    # println()
                    decreasing_step = false
                    accept = false
                end
            elseif max_displacement > solver_params.step_down_error # too large
                P.δt = P.δt / solver_params.stepping_factor
                decreasing_step = true
                # println()
                # println("Time step ↓ $(P.δt)")
                # println()
                if P.δt < solver_params.min_δt
                    throw("Time step became too small")
                end
                accept = false
            end

            if accept
                x .= x_candidate
                current_time += P.δt
                if F.centrosome
                    coords.centro_x .= centro_x_candidate
                    coords.centro_angle .= centro_angle_candidate
                    Centrosome.compute_vr(P, coords, centro_bufs, centro_vr)
                end
                break
            end
        end

        if metrics["started"]
            Metrics.update_metrics!(metrics, k, P, F, config, coords, coords_s, old_coords, nucleus_coords)
        end

        if F.cortex
            Cortex.copy(old_coords, coords)
        end
        if F.nucleus
            Nucleus.copy(old_nucleus_coords, nucleus_coords)
        end
        catch e
            metrics["crashed"] = true

            println()
            println("ERROR at iteration ", k, ":")
            println(e)
            break
        end

        # plot
        if (F.plot && ((k % plot_period == 0) || stepping))
            Plotting.update_plot(coords, nucleus_coords, k, current_time, P, F, false, plotables, centro_vr)
            if F.write_animation
                writer.grab_frame()
            end
            height = sum(x[:,2]/P.N)
            long_speed = (height - prev_height) / (plot_period*P.δt)
            prev_height = height
        end

        # writedlm("$(config["output_prefix"])Run_$(config["date_string"])_last_x.csv", coords.x, ',')
        # writedlm("$(config["output_prefix"])Run_$(config["date_string"])_last_centro_x.csv", coords.centro_x, ',')
    end

    Metrics.close_metrics!(metrics, k, current_time, P)
    Metrics.save_metrics(metrics, "$(config["output_prefix"])Run_$(config["date_string"])")

    if F.write_animation && F.plot
        writer.finish()
    end
    println("Finished")

    # if F.plot && metrics["started"]
        # Plotting.plot_metrics(metrics)
    # end

    if get(config, "batch", false)
        Plotting.close()
    end

    return metrics
end

end # module
