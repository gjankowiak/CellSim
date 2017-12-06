module Refactor

import AdhCommon
import Forces
import Wall
import Masks
import Utils
#import Surface
import SurfacesCommon
import EvenParam

import PyPlot
import NLsolve

function init_plot(x::Matrix, P::AdhCommon.Params, F::AdhCommon.Flags)
    PyPlot.ion()

    fig = PyPlot.figure()
    ax = PyPlot.axes(xlim = (-10, 10),ylim=(-15, 15))

    ax[:set_aspect]("equal", "datalim")

    line = ax[:plot](x[:,1], x[:,2], ".-", zorder=1)[1]
    #tracker = ax[:scatter](x[tr_idx,1], x[tr_idx,2], color="black", zorder=2)
    ax[:plot](x[:,1], x[:,2], color="black", lw=0.5)[1] # initial condition

    y = collect(linspace(-40, 40, 1000))
    wall = Wall.compute_walls(y, P)
    levelset = Wall.compute_walls(y, P, 1e-8)
    ax[:plot](wall, y, -wall, y, color="black", lw=0.5)

    levelset = Wall.compute_walls(y, P, 2e-4)
    ax[:plot](levelset, y, -levelset, y, color="red", lw=0.5)

    ax[:axvline](0)

    #quiver1 = ax[:quiver](x[:,1], x[:,2], quiver_data1[:,1], quiver_data1[:,2], units="x", width=0.001, color="blue")
    #quiver2 = ax[:quiver](x[:,1], x[:,2], quiver_data2[:,1], quiver_data2[:,2], scale_units="x", width=0.001, color="blue")
    #quiver3 = ax[:quiver](x[:,1], x[:,2], quiver_data3[:,1], quiver_data3[:,2], units="x", width=0.001, color="blue")
    #quiver4 = ax[:quiver](x[:,1], x[:,2], quiver_data4[:,1], quiver_data4[:,2], units="x", width=0.001, color="blue")

    #figManager = PyPlot.get_current_fig_manager()
    #figManager[:window][:showMaximized]()
    PyPlot.show()

    fig[:tight_layout]()
    return line
end

function update_plot(line, x::Matrix, k::Int, P::AdhCommon.Params, F::AdhCommon.Flags, initializing::Bool)
    ax = line[:axes]
    line[:set_data](x[:,1], x[:,2])
    #tracker[:set_data](x[tr_idx,1], x[tr_idx,2])

    #line_dbg[:set_data](1:N, x[:,1]-x0[:,1])

    # effective angle force
    #FA = @entry_norm(reshape(M_DFangle * δx, N, 2))/δt

    if initializing
        prefix = "[INIT]"
    else
        prefix = ""
    end
    ax[:set_title](@sprintf("%s N: %d, iter: %d", prefix, P.N, k))

    # effective angle force
    #FA = @entry_norm(reshape(M_DFangle * δx, N, 2))/δt

    #tracker[:set_offsets](x[tr_idx,:])
    #tracker[:set_sizes](3e3*FA[tr_idx])
    #tracker[:set_sizes](40*(drag_mask[tr_idx]/maximum(drag_mask)))
    #tracker[:set_sizes](tracker_data)
    #tracker[:set_facecolors](map(x-> x>1-1e-10 ? "#00FF00" : "#FF0000", drag_mask[tr_idx]))
    #tracker[:set_facecolors](map(x-> x>1-1e-10 ? "#00FF00" : "#FF0000", mask[tr_idx]))

    #quiver1[:set_offsets](x)
    #quiver1[:set_UVC](quiver_data1[:,1], quiver_data1[:,2])

    #quiver2[:set_offsets](x)
    #quiver2[:set_UVC](quiver_data2[:,1], quiver_data2[:,2])

    #quiver3[:set_offsets](x)
    #quiver3[:set_UVC](quiver_data3[:,1], quiver_data3[:,2])

    #quiver4[:set_offsets](x)
    #quiver4[:set_UVC](quiver_data4[:,1], quiver_data4[:,2])

    PyPlot.draw()
    sleep(0.0001)
    #read(STDIN, 1)
end

function main()
    # initialization

    N = 100

    AdhCommon.init(N)

    # Parameters
    P = AdhCommon.Params(
             N, # number of points
             100000, # max. number of iterations

             1/N, # space step
             1e-3, # time step

             3e-1, # pressure
             5e-2, # membrane elasticity
             1, # cortex viscosity

             2e-1, # polymerization speed

             16, # concentration of drag force in 1/(number of nodes)

             2*4.1, # initial ellipsis width
             2*4.1, # initial ellipsis height
             0.0, # initial vertical shift

             # Confinement field

             5e0, # sharpness
             0.0, # depth
             1, # pulsation
             1, # direction
             1, # number of Fourier components
             # to approximate a saw-tooth signal
             5.5, # mean width
             5.0 # inner width
            )

    println("equilibrium radius: ", 1/(2*pi - P.P/P.K))
    println("critical pressure: ", 2*pi*P.K)

    # Flags
    F = AdhCommon.Flags(
            false,  # confine
            false, # adjust_drag
            false, # polymerize
            true,  # dryrun
            true,  # plot
            false, # pretty
            true,  # continuous
            false, # innerloop
           )

    t = linspace(0, 1, N+1)[1:N]

    x_init = 0.5 * Float64[P.x0_a*cospi.(2t) P.x0_b*sinpi.(2t)]
    x = EvenParam.reparam(x_init, true)
    line = init_plot(x_init, P, F)

    Forces.init_FD_matrices(P)
    coords, coords_s = Forces.new_PointCoords(x, P)

    resi, resi_J = Forces.wrap_residuals(coords, coords_s, P, F)
    if F.innerloop
        resi_solver = NLsolve.DifferentiableSparseMultivariateFunction(resi, resi_J)
    else
        r_x = zeros(2N)
        Jr_x = spzeros(2N,2N)
        δx = zeros(2N)
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

        height = sum(x[:,2]/P.N)
        # println("Long. speed: ", (height - prev_height) / P.δt)
        prev_height = height

        # plot
        if F.plot & (k % 10 == 0)
            update_plot(line, x, k, P, F, false)
        end

        l2_norm = sqrt(sum(abs2, x))

        if l2_norm > 1e4
            println("Divergence detected, aborting")
            break
        end
    end

    println("Finished, type Enter to exit")
    read(STDIN, 1)
    # output
end

end # module
