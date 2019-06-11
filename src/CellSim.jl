module CellSim

import CellSimCommon
const CSC = CellSimCommon

import Cortex
import Nucleus
import Wall
import Masks

import Plotting

import Centrosome

import LinearAlgebra
const LA = LinearAlgebra

import JankoUtils
import EvenParam

import NLsolve

import YAML
import Dates

import DelimitedFiles: readdlm, writedlm

import SparseArrays
const SA = SparseArrays

macro eval_if_string(s)
    return esc(:(isa($s, String) ? eval(Meta.parse($s)) : $s))
end

function compute_initial_x(P::CSC.Params, F::CSC.Flags; convex::Bool=true)
    t = collect(range(0; stop=1, length=P.N+1))[1:P.N]

    if convex
        if !F.circular_wall
            return 0.5 * Float64[P.x0_a*cospi.(2t) P.x0_b*sinpi.(2t)]
        else
            return 0.5 * Float64[P.x0_a*cospi.(2t) .+ 2P.polar_shift P.x0_b*sinpi.(2t)]
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

    config_filename = "configs/default.yaml"
    if length(ARGS) > 0 && isfile(ARGS[1])
        config_filename = ARGS[1]
    else
        println("Please provide a configuration file!")
        println("Usage: run.jl <config>")
        exit(0)
    end

    date_string = string(Dates.now())

    yaml_config = YAML.load(open(config_filename))
    y_params = yaml_config["params"]

    output_prefix = "runs/"

    if haskey(yaml_config, "output_prefix")
        output_prefix = yaml_config["output_prefix"]
    end

    run(`mkdir -p runs`)
    run(`mkdir -p dumps`)
    run(`mkdir -p $(output_prefix)`)
    run(`cp $config_filename $(output_prefix)Run_$date_string.yaml`)


    # Parameters
    P = CSC.Params(
        y_params["M"],
        1/y_params["N"], # Δσ
        y_params["δt"],
        @eval_if_string(y_params["N"]),
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
        y_params["k_MT"],
        y_params["MT_potential_power"],
        y_params["MT_factor"],
        y_params["Nnuc"],
        y_params["N_P"],
        y_params["N_mu"],
        y_params["N_target_area"],
        y_params["N_kb"],
        y_params["N_ω"],
        y_params["N_W0"],
        y_params["N_kcont"],
        y_params["N_αcont"],
        y_params["N_kc"],
        y_params["N_l0c"],
        y_params["N_r_init"]
   )

    println("equilibrium radius: ", 1/(2*pi - P.P/P.K))
    println("critical pressure: ", 2*pi*P.K)

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
          y_flags["write_animation"],

          # debug
          y_flags["debug"]
    )

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

    if F.DEBUG
        println("*** DEBUG MODE ***")
    end

    run(`stty -icanon`)

    plotables = CSC.new_plotables(P.N)

    if haskey(yaml_config, "load_state")
        load_state = yaml_config["load_state"]
        if yaml_config["load_state"]["do_load"]
            x_init = readdlm(load_state["filename"], ',')
            if load_state["do_recenter"]
                x_init .-= sum(x_init; dims=1)/size(x_init, 1)
                #x_init[:,1] *= P.x0_a/abs(maximum(x_init[:,1]))
                #x_init[:,2] *= P.x0_b/abs(maximum(x_init[:,2]))
            end

            # if F.circular_wall
            # x_init[:,1] += P.polar_shift
            # end

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
            x_init = EvenParam.reparam(compute_initial_x(P, F; convex=true))
        end
    end


    x = copy(x_init)

    Cortex.init_FD_matrices(P)
    coords, coords_s = Cortex.new_PointCoords(x, P)

    if haskey(yaml_config, "load_state") && load_state["do_load"] && load_state["init_centro"]
        # if the initial centrosome location is given in the config, load it
        # the centrosome angle doesn't have any impact at this point
        coords.centro_x[:] = readdlm(load_state["filename_centro"], ',')
    else
        # otherwise, pick the centrosome location as the initial center of mass
        coords.centro_x[:] = sum(x_init; dims=1)/size(x_init,1)
    end

    nucleus_coords = missing

    if F.nucleus
        old_nucleus_coords = Nucleus.initialize_coords(P, F, coords)
        nucleus_coords = Nucleus.initialize_coords(P, F, coords)
        #
        # shift nucleus to a wide section
        #
        if haskey(yaml_config, "load_state") && haskey(load_state, "shift_nucleus")
            old_nucleus_coords.Y[:] = old_nucleus_coords.Y .- sum(old_nucleus_coords.Y; dims=1)/P.Nnuc .+ [0.0 load_state["shift_nucleus"]]
            nucleus_coords.Y[:] = nucleus_coords.Y .- sum(nucleus_coords.Y; dims=1)/P.Nnuc .+ [0.0 load_state["shift_nucleus"]]
        end

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

    # centrosome buffers and coordinates
    (centro_bufs, centro_vr, centro_qw, centro_pc) = Centrosome.init(P)
    Centrosome.compute_vr(P, coords, centro_bufs, centro_vr)
    (centro_A, centro_id_comp, centro_b_ce, centro_b_ce_rhs, centro_b_co_rhs) = Centrosome.assemble_system(P, F, coords, centro_bufs, centro_vr, centro_qw, centro_pc, plotables, potentials)

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
        fig = Plotting.init_plot(coords, P, F)
        Plotting.update_plot(coords, nucleus_coords, 0, P, F, false, plotables, centro_vr)
        if F.write_animation
            writer = Plotting.init_animation(date_string)
            writer.setup(fig, string(writer.metadata["title"], ".mp4"), 100)
        end
    end

    k = 0
    prev_height = 0.0

    metrics = Dict{String, Float64}()

    stepping = F.DEBUG
    if !stepping
        input_task = @async read(stdin, Char)
    end

    println("Initialized")

    # outer loop
    while k < P.M

        if stepping
            println("debug: 'q' to quit, 'c' to run continuously, any other key to step, 'b' to run step by step")
            key = read(stdin, 1)[1]
            if key == 0x63 # c
                stepping = false
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

        k += 1
        print("\b"^100)
        println(" iteration #", k, ", ")

        # initialization metrics
        if k == 400
            metrics["t_start"] = k*P.δt
            barycenter_y = sum(x[:,2])/P.N
            metrics["barycenter_y_start"] = barycenter_y
        end


        # inner loop
        if k > 1
            Cortex.update_coords(coords, P, x)
        end

        if F.nucleus
            fill!(potentials.N_W, 0.0)
            fill!(potentials.N_∇W, 0.0)
            fill!(potentials.C_∇W, 0.0)
            fill!(potentials.CS_∇W, 0.0)
            Nucleus.compute_contact_force(potentials, coords, nucleus_coords, P, F)
            if F.centrosome
                Nucleus.compute_centronuclear_force(potentials, coords, nucleus_coords, P, F)
            end
            Nucleus.update_coords(old_nucleus_coords, nucleus_coords, potentials, P, F, temparrays)

            Nucleus.copy(old_nucleus_coords, nucleus_coords)
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

            x .+= reshape(δx[1:2P.N], P.N, 2)

            coords.centro_x .+= δx[(2P.N+1):(2P.N+2)]
            coords.centro_angle .+= δx[2P.N+3]
            Centrosome.compute_vr(P, coords, centro_bufs, centro_vr)
        else
            δx[:] = -Jr_x\r_x
            x .+= + reshape(δx, P.N, 2)
        end

        # plot
        plot_period = (F.write_animation || F.DEBUG) ? 1 : 1
        if (F.plot & (k % plot_period == 0))
            Plotting.update_plot(coords, nucleus_coords, k, P, F, false, plotables, centro_vr)
            if F.write_animation
                writer.grab_frame()
            end
            height = sum(x[:,2]/P.N)
            long_speed = (height - prev_height) / (plot_period*P.δt)
            prev_height = height
            # println("Long. speed: ", long_speed)
        end

        writedlm("$(output_prefix)Run_$(date_string)_last_x.csv", coords.x, ',')
        writedlm("$(output_prefix)Run_$(date_string)_last_centro_x.csv", coords.centro_x, ',')

        l2_norm = sqrt(sum(abs2, x))

        # if l2_norm > 1e4
        # println("Divergence detected, aborting")
        # break
        # end


    end

    metrics["t_end"] = k*P.δt
    barycenter_y = sum(x[:,2])/P.N
    metrics["barycenter_y_end"] = barycenter_y
    metrics["speed"] = (metrics["barycenter_y_end"] - metrics["barycenter_y_start"])/(metrics["t_end"] - metrics["t_start"])
    metrics["fw0"] = P.f_ω0
    metrics["fb"] = P.f_β
    metrics["fwidth"] = P.f_width
    JankoUtils.write_dict("$(output_prefix)Run_$(date_string)_metrics.yaml", metrics)

    if F.write_animation & F.plot
        writer.finish()
    end
    println("Finished")

    return
end

end # module
