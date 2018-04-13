const sqrtEPS = 1e-5

import Forces: PointCoords, PointCoordsShifted
import CellSimCommon

function compute_angles_radii(vr::Array{Float64,2})
    n = size(vr, 1)

    # FIXME
    idx_p = circshift(1:n, -1)
    θ = angles.(vr[:,1] + vr[:,2]*im)
    θ_p = view(θ, idx_p)

    r = CellSimCommon.@entry_norm(vr)
    r_p = view(r, idx_p)

    return θ, θ_p, r, r_p
end

function compute_poly_coefficients(θ::Vector{Float64}, θ_p::Vector{Float64}, r::Vector{Float64}, r_p)
    a = @. (r_p*sin(θ_p) - r*sin(θ))/(r*r_p*sin(θ_p-θ))
    b = @. (-r_p*cos(θ_p) - r*cos(θ))/(r*r_p*sin(θ_p-θ))
    return a, b
end

function quadrature_weights(coords::PointCoords, coords_s::PointCoordsShifted, vr::Array{Float64,2}, Δvr::Array{Float64,2}, z::Vector{Float64})
    c0 = sum(abs2, vr .- reshape(z, 1, 2), 2)
    c1 = 2*CellSimCommon.@dotprod(coords.Δvr, coords.vr - reshape(z, 1, 2))
    c2 = sum(abs2, coords.ΔL.^2, 2)

    delta = sqrt.(4*c2.*c0-c1.^2)

    I0 = (delta .> sqrtEPS) .* (2*(atan.((2c2+c1)./delta)-atan.(c1./delta))./delta)
    I1 = (-c1.*I0 + log.(1+(c1+c2)./c0))./2c2
    I2 = -c1./c2*I1-c0./c2.*I0+1./c2

    w1 = Δσ*
