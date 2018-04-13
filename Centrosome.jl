module Centrosome

const sqrtEPS = 1e-5

import Forces: PointCoords, PointCoordsShifted
import CellSimCommon
import Visibility

mutable struct VisibleRegion
    nodes::Array{Float64,2}
    diffs::Array{Float64,2}
    M_inter::SparseMatrixCSC
    Δs::Vector{Float64}
    n::Int64
end

struct QuadratureWeights
    c0::Vector{Float64}
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
    θ_p::SubArray{Float64,1}
    r_p::SubArray{Float64,1}
end

function compute_angles_radii(vr::Array{Float64,2}, polar::PolarCoordinates)
    n = size(vr.n, 1)

    # FIXME
    idx_p = circshift(1:n, -1)
    polar.θ[:] = angles.(vr[:,1] + vr[:,2]*im)
    polar.θ_p = view(polar.θ, idx_p)

    polar.r[:] = CellSimCommon.@entry_norm(vr)
    polar.r_p = view(polar.r, idx_p)
end

function compute_line_coefficients(polar::PolarCoordinates, qw::QuadratureWeights)
    qw.a[:] = @. (polar.r_p*sin(polar.θ_p) - polar.r*sin(polar.θ))/(polar.r*polar.r_p*sin(polar.θ_p-polar.θ))
    qw.b[:] = @. (-polar.r_p*cos(polar.θ_p) - polar.r*cos(polar.θ))/(polar.r*polar.r_p*sin(polar.θ_p-polar.θ))
end

function quadrature_weights(vr::VisibleRegion, qw::QuadratureWeights)
    idx_s = circshift(1:vr.n, 1)

    qw.c0[:] = sum(abs2, vr.nodes, 2)
    qw.c1[:] = 2*CellSimCommon.@dotprod(vr.diffs, vr.nodes)
    qw.c2[:] = sum(abs2, vr.diffs, 2)

    qw.delta[:] = sqrt.(4*qw.c2.*qw.c0-qw.c1.^2)

    qw.I0[:] = (qw.delta .> sqrtEPS) .* (2*(atan.((2qw.c2+qw.c1)./qw.delta)-atan.(qw.c1./qw.delta))./qw.delta)
    qw.I1[:] = -qw.c1.*qw.I0 + log.(1+(qw.c1+qw.c2)./qw.c0))./2qw.c2
    qw.I2[:] = -qw.c1./qw.c2*qw.I1-qw.c0./qw.c2.*qw.I0+1./qw.c2

    @views qw.w2[:] = abs.(-vr.diffs[:,2].*vr.nodes[:,1]+vr.diffs[:,1].*vr.nodes[:,2])
    qw.w1[:] = vr.Δs.*qw.w2.*(qw.I0-qw.I1)

    qw.wp1[:] = vr.Δs.*qw.w2.*(qw.I1-qw.I2)
    qw.wp2[:] = vr.Δs.*qw.w2.*qw.I2

    qw.w2[:] .*= vr.Δs.*qw.I1

    qw.A1[:] = qw.w1 + qw.w2[idx_s]

    @views qw.A2[1,:] = @. - qw.w1*vr.nodes[:,2] - qw.wp1*vr.diffs[:,2] - qw.w2[idx_s]*vr.nodes[idx_s,2] - qw.wp2[idx_s]*vr.diffs[idx_s,2]
    @views qw.A2[2,:] = @. qw.w1*vr.nodes[:,1] + qw.wp1*vr.diffs[:,1] + qw.w2[idx_s]*vr.nodes[idx_s,1] + qw.wp2[idx_s]*vr.diffs[idx_s,1]
end

function compute_vr(P::Params, coords::PointCoords, bufs::Visibility.RawVisibleRegion, vr::VisibleRegion)
    raw_vr = Visibility.vispol(coords.x, coords.centro, bufs)
    vr.nodes[1:raw_vr.t,:] = raw_vr.q[1:raw_vr.t,:]
    circshift!(vr.diffs, vr.nodes, -1)
    vr.diffs += vr.nodes
    vr.n = raw_vr.t
    vr.M_inter = Visibility.build_interpolation_matrix(raw_vr, P.N)
    vr.Δs = vr.M_inter * (P.Δσ*(1:P.N))
    @views vr.Δs[1:vr.n] = -vr.Δs[1:vr.n] + circshift(vr.Δs[1:vr.n], -1)
end

# FIXME add initialization

end # module Centrosome
