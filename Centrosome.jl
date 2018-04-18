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
                           zeros(P.N+1), zeros(P.N+1), zeros(2, P.N+1), # delta, A1, A2
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
    @views qw.A2[1,1:vr.n] = @. - qw.w1[1:vr.n]*vr.nodes[1:vr.n,2] - qw.wp1[1:vr.n]*vr.diffs[1:vr.n,2] - qw.w2[idx_m]*vr.nodes[idx_m,2] - qw.wp2[idx_m]*vr.diffs[idx_m,2]
    @views qw.A2[2,1:vr.n] = @. qw.w1[1:vr.n]*vr.nodes[1:vr.n,1] + qw.wp1[1:vr.n]*vr.diffs[1:vr.n,1] + qw.w2[idx_m]*vr.nodes[idx_m,1] + qw.wp2[idx_m]*vr.diffs[idx_m,1]
    qw.A2[(vr.n+1):end,:] = 0.0
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
    p = 0.0
    k = 5e0
    plotables.mt_force_indiv[:] = k*vr.nodes.*sign.(CellSimCommon.@entry_norm(vr.nodes)-MT_RADIUS).*abs.(CellSimCommon.@entry_norm(vr.nodes)-MT_RADIUS).^p

    # pulling
    # p = 0.0
    # k = 5e0
    # plotables.mt_force_indiv[:] = k*vr.nodes.*CellSimCommon.@entry_norm(vr.nodes).^p

    # pushing
    # p = -1
    # k = 5e0
    # plotables.mt_force_indiv[:] = -k*vr.nodes.*abs.(CellSimCommon.@entry_norm(vr.nodes)-MT_RADIUS).^p
end

function integrate_θ(f::Array{Float64}, qw::QuadratureWeights, vr::VisibleRegion)
    return sum(f[1:vr.n,:].*qw.A1[1:vr.n], 1)
end

function interpolate_on_nodes(f::Array{Float64}, vr::VisibleRegion)
    return vr.M_inter*f
end

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
    @views A[3,3] = CellSimCommon.nansum(@. (pc.θ_p[1:n] - pc.θ[1:n])*(
                                                1/(qw.a[1:n]^2 * cot(pc.θ[vr.idx_p]) + qw.a[1:n]*qw.b[1:n])
                                                -1/(qw.a[1:n]^2 * cot(pc.θ[1:n]) + qw.a[1:n]*qw.b[1:n])))

    @views A[1,3] = CellSimCommon.nansum(@. (pc.θ_p[1:n] - pc.θ[1:n])/(qw.a[1:n]^2 + qw.b[1:n]^2)*(
                                            (qw.a[1:n]*log(pc.r[vr.idx_p]) - qw.b[1:n]*pc.θ[vr.idx_p])
                                            -(qw.a[1:n]*log(pc.r[1:n]) + qw.b[1:n]*pc.θ[1:n])))

    @views A[2,3] = CellSimCommon.nansum(@. (pc.θ_p[1:n] - pc.θ[1:n])/(qw.a[1:n]^2 + qw.b[1:n]^2)*(
                                            (qw.b[1:n]*log(pc.r[vr.idx_p]) + qw.a[1:n]*pc.θ[vr.idx_p])
                                            -(qw.b[1:n]*log(pc.r[1:n]) - qw.a[1:n]*pc.θ[1:n])))
    A[3,1] = A[1,3]
    A[3,2] = A[2,3]
    A *= -P.k_MT

    b = zeros(3)

    compute_mt_longi_force(vr, plotables)

    b[1:2] = integrate_θ(plotables.mt_force_indiv, qw, vr)

    return A, b
end

end # module Centrosome
