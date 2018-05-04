module Centrosome

const sqrtEPS = 1e-5
const EPS = sqrtEPS^2

const MT_RADIUS = 0.05

import Forces: PointCoords, PointCoordsShifted
import CellSimCommon
import Visibility

#=
Except indices vectors, all vectors here are preallocated
Since the visibility region doesn't have a constant number
of nodes, this means that part of the end of the vectors
are meaningless.

The meaningful part is between indices 1 and vr.n, after the
visibility region has been computed with "compute_vr"
=#

mutable struct VisibleRegion
    nodes::Array{Float64,2}   # nodes of the visibility region
                              # RELATIVE TO THE CENTROSOME
    diffs::Array{Float64,2}   # nodes[i+1] - nodes[i]
    lengths::Vector{Float64}  # |diffs[i]|
    M_inter::SparseMatrixCSC  # interpolation matrix f(X) -> f(nodes)
    Δs::Vector{Float64}       # s(i+1) - s(i)
    n::Int64                  # number of nodes in the visibility region
    idx_m::Vector{Int64}      # index of the preview node (_m for -, minus)
    idx_p::Vector{Int64}      # index of the next node (_p for +, plus)
    lambdas::Vector{Float64}  # (1-λ) where λ is the convex combination parameter on each edge
end

struct QuadratureWeights
    c0::Vector{Float64}       # see articles appendix
    c1::Vector{Float64}
    c2::Vector{Float64}
    I0::Vector{Float64}
    I1::Vector{Float64}
    I2::Vector{Float64}
    w1::Vector{Float64}
    w2::Vector{Float64}
    wp1::Vector{Float64}
    wp2::Vector{Float64}
    delta::Vector{Float64}
    A1::Vector{Float64}
    A2::Array{Float64,2}
    a::Vector{Float64}
    b::Vector{Float64}
end

mutable struct PolarCoordinates
    θ::Vector{Float64}
    r::Vector{Float64}
    θ_p::Vector{Float64}
    r_p::Vector{Float64}
end

function init(P::CellSimCommon.Params)
    bufs = Visibility.init(P.N+1)
    vr = VisibleRegion(zeros(P.N+1, 2), zeros(P.N+1, 2), zeros(P.N+1),
                       spzeros(0,0), zeros(P.N+1), 0, Int64[], Int64[], zeros(P.N+1))
    qw = QuadratureWeights(zeros(P.N+1), zeros(P.N+1), zeros(P.N+1), # c0, c1, c2
                           zeros(P.N+1), zeros(P.N+1), zeros(P.N+1), # I0, I1, I2,
                           zeros(P.N+1), zeros(P.N+1), # w1, w2
                           zeros(P.N+1), zeros(P.N+1), # wp1, wp2
                           zeros(P.N+1), zeros(P.N+1), zeros(P.N+1, 2), # delta, A1, A2
                           zeros(P.N+1), zeros(P.N+1)) # a, b
    θ = zeros(P.N+1)
    r = zeros(P.N+1)
    pc = PolarCoordinates(θ, r, view(θ, :), view(r, :))

    return (bufs, vr, qw, pc)
end

function compute_angles_radii(vr::VisibleRegion, polar::PolarCoordinates)
    polar.θ[:] = angle.(vr.nodes[:,1] + vr.nodes[:,2]*im)
    polar.θ_p[1:vr.n] = polar.θ[vr.idx_p]

    polar.r[:] = CellSimCommon.@entry_norm(vr.nodes)
    polar.r_p[1:vr.n] = polar.r[vr.idx_p]
end

function compute_line_coefficients(polar::PolarCoordinates, qw::QuadratureWeights)
    qw.a[:] = @. (polar.r_p*sin(polar.θ_p) - polar.r*sin(polar.θ))/(polar.r*polar.r_p*sin(polar.θ_p-polar.θ))
    qw.b[:] = @. (-polar.r_p*cos(polar.θ_p) - polar.r*cos(polar.θ))/(polar.r*polar.r_p*sin(polar.θ_p-polar.θ))
end

function compute_quadrature_weights(vr::VisibleRegion, qw::QuadratureWeights)
    # FIXME restrict computations to 1:vr.n
    # doing it on the whole range is correct but results some useless computations
    # this is fine for a "small" blind region

    qw.c0[:] = sum(abs2, vr.nodes, 2)
    qw.c1[:] = 2*CellSimCommon.@dotprod(vr.diffs, vr.nodes)
    qw.c2[:] = sum(abs2, vr.diffs, 2)

    qw.delta[:] = sqrt.(max.(4*qw.c2.*qw.c0-qw.c1.^2, 0.0))
    idx = find(qw.delta .> sqrtEPS)

    fill!(qw.I0, 0.0)
    fill!(qw.I1, 0.0)
    fill!(qw.I2, 0.0)

    @views begin
        qw.I0[idx] = (2*(atan.((2qw.c2[idx] + qw.c1[idx])./qw.delta[idx]) - atan.(qw.c1[idx]./qw.delta[idx]))./qw.delta[idx])
        qw.I1[idx] = @. (-qw.c1[idx]*qw.I0[idx]  +  log(1 + (qw.c1[idx] + qw.c2[idx])/qw.c0[idx]))/2qw.c2[idx]
        qw.I2[idx] = @. -qw.c1[idx]/qw.c2[idx]*qw.I1[idx] - qw.c0[idx]/qw.c2[idx]*qw.I0[idx] + 1/qw.c2[idx]
    end

    @views qw.w2[:] = abs.(-vr.diffs[:,2].*vr.nodes[:,1]+vr.diffs[:,1].*vr.nodes[:,2])

    qw.w1[:] = vr.Δs.*qw.w2.*(qw.I0-qw.I1)

    qw.wp1[:] = vr.Δs.*qw.w2.*(qw.I1-qw.I2)
    qw.wp2[:] = vr.Δs.*qw.w2.*qw.I2

    qw.w2[:] .*= vr.Δs.*qw.I1

    qw.A1[(vr.n+1):end] = 0.0
    qw.A1[1:vr.n] = qw.w1[1:vr.n] + qw.w2[vr.idx_m]

    idx_m = vr.idx_m
    @views qw.A2[1:vr.n,1] = @. - qw.w1[1:vr.n]*vr.nodes[1:vr.n,2] - qw.wp1[1:vr.n]*vr.diffs[1:vr.n,2] - qw.w2[idx_m]*vr.nodes[idx_m,2] - qw.wp2[idx_m]*vr.diffs[idx_m,2]
    @views qw.A2[1:vr.n,2] = @. qw.w1[1:vr.n]*vr.nodes[1:vr.n,1] + qw.wp1[1:vr.n]*vr.diffs[1:vr.n,1] + qw.w2[idx_m]*vr.nodes[idx_m,1] + qw.wp2[idx_m]*vr.diffs[idx_m,1]
    qw.A2[:,(vr.n+1):end] = 0.0
end

function compute_vr(P::CellSimCommon.Params, coords::PointCoords, bufs::Visibility.Buffers, vr::VisibleRegion)
    raw_vr = Visibility.vispol(coords.x, coords.centro_x, bufs)

    vr.n = raw_vr.t

    idx_m = circshift(1:vr.n, 1)
    vr.idx_m = idx_m
    idx_p = circshift(1:vr.n, -1)
    vr.idx_p = idx_p

    # the visibility region nodes are given RELATIVE TO THE CENTROSOME
    vr.nodes[1:vr.n,:] = raw_vr.q[1:vr.n,:] .- reshape(coords.centro_x, 1, 2)
    @views circshift!(vr.diffs[1:vr.n,:], vr.nodes[1:vr.n,:], -1)
    vr.diffs += vr.nodes
    vr.lengths = CellSimCommon.@entry_norm(vr.diffs)

    vr.lambdas[1:vr.n] = raw_vr.lambdas[1:vr.n]

    vr.M_inter = Visibility.build_interpolation_matrix(raw_vr, P.N)

    vr.Δs[1:vr.n] = P.Δσ * vr.lambdas[1:vr.n]
end

function compute_mt_longi_force(vr::VisibleRegion, plotables::CellSimCommon.Plotables)

    # pushing upto MT_RADIUS and the pulling
    # p = 1.0
    # k = 5e0
    # plotables.mt_force_indiv[:] = k*vr.nodes./CellSimCommon.@entry_norm(vr.nodes).*sign.(CellSimCommon.@entry_norm(vr.nodes)-MT_RADIUS).*abs.(CellSimCommon.@entry_norm(vr.nodes)-MT_RADIUS).^p

    # pulling
    #p = 0.0
    #k = 5e0
    #plotables.mt_force_indiv[:] = k*vr.nodes.*CellSimCommon.@entry_norm(vr.nodes).^p

    # pushing
    p = -2
    k = 1e0
    plotables.mt_force_indiv[1:vr.n,:] = -k*vr.nodes[1:vr.n,:].*abs.(CellSimCommon.@entry_norm(vr.nodes[1:vr.n,:])).^p
end

function integrate_θ(f::Array{Float64}, qw::QuadratureWeights, vr::VisibleRegion)
    return sum(f[1:vr.n,:].*qw.A1[1:vr.n], 1)
end

function integrate_θ_eθ(f::Array{Float64}, qw::QuadratureWeights, vr::VisibleRegion)
    return sum(qw.A2[1:vr.n,:].*f[1:vr.n,:])
end

function interpolate_on_nodes(f::Array{Float64}, vr::VisibleRegion)
    return vr.M_inter*f
end

"""
The full system (cortex + centrosome) is of the form

[[ M    | b_co ]    (X_co^(n+1))    (b_co_n X_ce^n)
 [ b_ce | A    ]]   (X_ce^(n+1)) =  (b_ce_n X_co^n)
"""
function assemble_system(P::CellSimCommon.Params, coords::PointCoords, bufs::Visibility.Buffers, vr::VisibleRegion,
                        qw::QuadratureWeights, pc::PolarCoordinates, plotables::CellSimCommon.Plotables)
    compute_vr(P, coords, bufs, vr)
    compute_angles_radii(vr, pc)
    compute_line_coefficients(pc, qw)
    compute_quadrature_weights(vr, qw)

    n = vr.n

    # for now we assume ∂t X = 0
    A = 2π*eye(3)

    # maybe wrap this into a function
    #= compute using analytical formula
    @views A[3,3] = CellSimCommon.nansum(@. ( 1/(qw.a[1:n]^2 * cot(pc.θ[vr.idx_p]) + qw.a[1:n]*qw.b[1:n])
                                             -1/(qw.a[1:n]^2 * cot(pc.θ[1:n]) + qw.a[1:n]*qw.b[1:n])))

    @views A[1,3] = CellSimCommon.nansum(@. 1/(qw.a[1:n]^2 + qw.b[1:n]^2)*(
                                             (qw.a[1:n]*log(pc.r[vr.idx_p]) - qw.b[1:n]*pc.θ[vr.idx_p])
                                            -(qw.a[1:n]*log(pc.r[1:n])      + qw.b[1:n]*pc.θ[1:n])))

    @views A[2,3] = CellSimCommon.nansum(@. 1/(qw.a[1:n]^2 + qw.b[1:n]^2)*(
                                             (qw.b[1:n]*log(pc.r[vr.idx_p]) + qw.a[1:n]*pc.θ[vr.idx_p])
                                            -(qw.b[1:n]*log(pc.r[1:n])      - qw.a[1:n]*pc.θ[1:n])))
    =#

    #= compute using quadrature =#
    @views A[3,3] = integrate_θ(pc.r[1:n].^2, qw, vr)[1]
    @views A[1:2,3] = integrate_θ([-vr.nodes[1:n,2] vr.nodes[1:n,1]], qw, vr)

    A[3,1] = A[1,3]
    A[3,2] = A[2,3]
    A *= P.k_MT

    b_ce = zeros(3, 2P.N)
    b_ce_n = zeros(3)

    # -k_MT/δt ∫ Xθ^(n+1)
    b_ce[1,1:P.N] = b_ce[2,(P.N+1):2P.N] = -P.k_MT*reshape(qw.A1[1:vr.n], 1, vr.n)*vr.M_inter

    # -k_MT/δt ∫ |Xc - Xθ| eθ^T . Xθ^(n+1)
    b_ce[3,1:P.N] = -P.k_MT*reshape(qw.A2[1:vr.n,1], 1, vr.n)*vr.M_inter
    b_ce[3,(P.N+1):2P.N] = -P.k_MT*reshape(qw.A2[1:vr.n,2], 1, vr.n)*vr.M_inter

    compute_mt_longi_force(vr, plotables)

    # ∫ F_MT
    b_ce_n[1:2] = integrate_θ(plotables.mt_force_indiv, qw, vr)

    #  - k_MT ∫ ds/dt ∂_s Xθ^n
    b_ce_n[1:2] +=  - P.k_MT*reshape(integrate_θ(vr.M_inter*reshape(plotables.transport_force, P.N, 2), qw, vr), 2, 1)

    # - k_MT ∫ |Xc - Xθ| eθ^T . ds/dt ∂_s Xθ^n
    b_ce_n[3] = -P.k_MT*integrate_θ_eθ(vr.M_inter*reshape(plotables.transport_force, P.N, 2), qw, vr)


    ## b_co = b_ce'
    # b_co = zeros(2P.N, 3)

    b_co_n = zeros(2P.N)
    b_co_n[1:P.N] = -vr.M_inter'*(qw.A1[1:vr.n].*plotables.mt_force_indiv[1:vr.n,1])
    b_co_n[(1+P.N):2P.N] = -vr.M_inter'*(qw.A1[1:vr.n].*plotables.mt_force_indiv[1:vr.n,2])

    b_co_n[1:P.N] += -P.k_MT*(vr.M_inter'*(qw.A1[1:vr.n].*(vr.M_inter*plotables.transport_force[1:P.N])))
    b_co_n[(P.N+1):2P.N] += -P.k_MT*(vr.M_inter'*(qw.A1[1:vr.n].*(vr.M_inter*plotables.transport_force[(P.N+1):2P.N])))

    id_comp = -P.k_MT/P.δt * spdiagm(repmat(vr.M_inter'*(qw.A1[1:vr.n].*(vr.M_inter*ones(P.N))), 2))

    # DEBUG, needs calibration
    #plotables.mt_force_indiv[1:vr.n,:] .*= 1e3*qw.A1[1:vr.n,:]

    return A, id_comp, b_ce, b_ce_n, b_co_n
end

end # module Centrosome
