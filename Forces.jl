module Forces

using AdhCommon
import Wall, Masks
import Utils: spdiagm_const, idxmod

immutable PointCoords
    x::Matrix{Float64}
    Δx::Matrix{Float64}
    Δ2x::Matrix{Float64}
    Δ2x_perp::Matrix{Float64}

    ΔL::Vector{Float64}
    Δ2L::Vector{Float64}

    ell::Vector{Float64}

    # unit tangent
    τ::Matrix{Float64}

    # centered unit tangent
    τc::Matrix{Float64}

    # unit normal
    v::Matrix{Float64}

    # centered unit tangent
    vc::Matrix{Float64}

    # derivatives of tangent and normal unit vectors
    dτc::Matrix{Float64}
    d2τ::Matrix{Float64}
    d3τ::Matrix{Float64}

    dvc::Vector{Float64}

    # auxiliary quantities
    dvΔL::Vector{Float64}
    vd2τ::Vector{Float64}
    d1vd2τ::Vector{Float64}
    d2vd2τ::Vector{Float64}
end

immutable PointCoordsShifted
    Δx_m::SubArray
    Δx_p::SubArray
    ΔL_m::SubArray
    ΔL_p::SubArray
    ell_m::SubArray
    τ_m::SubArray
    v_m::SubArray
    v_p::SubArray
end

type Differentials
    Dτ::SparseMatrixCSC{Float64}
    Dτc::SparseMatrixCSC{Float64}
    Dv::SparseMatrixCSC{Float64}
    Dvc::SparseMatrixCSC{Float64}
    Dτ_m::SparseMatrixCSC{Float64}
end

function new_PointCoords(N::Int64)
    coords = PointCoords(
         zeros(N,2),  # x
         zeros(N,2),  # Δx
         zeros(N,2),  # Δ2x
         zeros(N,2),  # Δ2x_perp
         zeros(N),  # ΔL
         zeros(N),  # Δ2L
         zeros(N),  # ell
         zeros(N,2),  # τ
         zeros(N,2),  # τc
         zeros(N,2),  # v
         zeros(N,2),  # vc
         zeros(N,2),  # dτc
         zeros(N,2),  # d2τ
         zeros(N,2),  # d3τ
         zeros(2N), # dvc
         zeros(2N), # dvΔL
         zeros(N), # vd2τ
         zeros(N), # d1vd2τ
         zeros(N)) # d2vd2τ

    coords_shifted = new_PointCoordsShifted(coords)
    return coords, coords_shifted
end

function new_PointCoords(x::Matrix{Float64}, P::Params)
    N = size(x, 1)
    coords, coords_shifted = new_PointCoords(N)
    update_coords(coords, P, x)
    return coords, coords_shifted
end

function new_PointCoordsShifted(coords::PointCoords)
    N = size(coords.x, 1)

    return PointCoordsShifted(
        view(coords.Δx, AdhCommon.circ_idx_m1, :),
        view(coords.Δx, AdhCommon.circ_idx_p1, :),
        view(coords.ΔL, AdhCommon.circ_idx_m1),
        view(coords.ΔL, AdhCommon.circ_idx_p1),
        view(coords.ell, AdhCommon.circ_idx_m1),
        view(coords.τ, AdhCommon.circ_idx_m1, :),
        view(coords.v, AdhCommon.circ_idx_m1, :),
        view(coords.v, AdhCommon.circ_idx_p1, :))
end

function new_Differentials(coords::PointCoords, coords_s::PointCoordsShifted)
    Dτ = AdhCommon.@bc_scalar(1./coords.ΔL).*(AdhCommon.pointwise_projection(coords.v))*D1p_unorm
    Dτc = AdhCommon.@bc_scalar(1./coords.Δ2L).*AdhCommon.pointwise_projection(coords.vc)*D1c_unorm
    Dv = -AdhCommon.@bc_scalar(1./coords.ΔL).*(AdhCommon.pointwise_projection(coords.τ))*D1p_perp_unorm
    Dvc = -AdhCommon.@bc_scalar(1./coords.Δ2L).*AdhCommon.pointwise_projection(coords.τc)*D1c_perp_unorm
    return Differentials(
        Dτ, Dτc, Dv, Dvc, M_cs_minus*Dτ
    )
end

function update_coords(coords::PointCoords, P::Params, x::Matrix{Float64})
    coords.x[:] = x
    AdhCommon.@delta!(x, coords.Δx)
    coords.Δ2x[:] = x[AdhCommon.circ_idx_p1,:] - x[AdhCommon.circ_idx_m1,:]
    coords.Δ2x_perp[:,1] = -coords.Δ2x[:,2]
    coords.Δ2x_perp[:,2] =  coords.Δ2x[:,1]

    coords.ΔL[:] = AdhCommon.@entry_norm(coords.Δx)
    coords.Δ2L[:] = AdhCommon.@entry_norm(coords.Δ2x)

    coords.ell[:] = coords.ΔL/P.Δσ - 1

    # unit tangent
    coords.τ[:] = coords.Δx./coords.ΔL

    # centered unit tangent
    coords.τc[:] = coords.Δ2x./AdhCommon.@entry_norm(coords.Δ2x)

    # unit normal
    coords.v[:] = hcat(view(coords.τ,:,2), -view(coords.τ,:,1))

    # centered unit tangent
    coords.vc[:] = -coords.Δ2x_perp./AdhCommon.@entry_norm(coords.Δ2x)

    # derivatives of tangent and normal unit vectors
    coords.dτc[:] = D1c*vec(coords.τc)
    coords.d2τ[:] = D2*vec(coords.τc)
    coords.d3τ[:] = D3*vec(coords.τc)

    coords.dvc[:] = D1c*vec(coords.vc)

    # auxiliary quantities
    coords.dvΔL[:] = 2P.Δσ*D1c_short*(coords.vc./coords.Δ2L)
    coords.vd2τ[:] = AdhCommon.@dotprod(coords.vc, coords.d2τ)
    coords.d1vd2τ[:] = D1c_short*coords.vd2τ
    coords.d2vd2τ[:] = D2_short*coords.vd2τ
end

function init_FD_matrices(P::Params)
    N = P.N
    Δσ = P.Δσ

    const global M_perp = blkdiag(-speye(N), speye(N))
    const global M_cs_plus = spdiagm(([1; zeros(N-1); 1], [ones(N-1); 0; ones(N-1)]), (-N+1, 1), 2N, 2N)
    const global M_cs_minus = spdiagm(([1; zeros(N-1); 1], [ones(N-1); 0; ones(N-1)]), (N-1, -1), 2N, 2N)

    const global D1p_unorm = blkdiag(spdiagm((1, -ones(N), ones(N-1)), (-N+1, 0, 1), N, N),
                        spdiagm((1, -ones(N), ones(N-1)), (-N+1, 0, 1), N, N))
    const global D1p = D1p_unorm/Δσ
    const global D1m = blkdiag(spdiagm((-ones(N-1), ones(N), -1), (-1, 0, N-1), N, N),
                         spdiagm((-ones(N-1), ones(N), -1), (-1, 0, N-1), N, N))/Δσ

    const global D1p_perp_unorm = ([[spzeros(N,N) -spdiagm((1, -ones(N), ones(N-1)), (-N+1, 0, 1), N, N)]
                       [spdiagm((1, -ones(N), ones(N-1)), (-N+1, 0, 1), N, N) spzeros(N,N)]])
    const global D1p_perp = D1p_perp_unorm/Δσ

    const global D1m_perp = ([[spzeros(N,N) -spdiagm((-ones(N-1), ones(N), -1), (-1, 0, N-1), N, N)]
                       [spdiagm((-ones(N-1), ones(N), -1), (-1, 0, N-1), N, N) spzeros(N,N)]])/Δσ

    # D1 matrix, centered differences
    # for point centered normal and tangent
    const global D1c_short_unorm = spdiagm((1, -ones(N-1), ones(N-1), -1), (-N+1, -1, 1, N-1), N, N)
    const global D1c_unorm       = blkdiag(D1c_short_unorm, D1c_short_unorm)
    const global D1c_perp_unorm  = ([[spzeros(N,N) -D1c_short_unorm];[D1c_short_unorm spzeros(N,N)]])
    const global D1c_short       = D1c_short_unorm/(2Δσ)
    const global D1c             = D1c_unorm/(2Δσ)
    const global D1c_perp        = D1c_perp_unorm/(2Δσ)

    # D2 matrix, centered differences
    # for point centered normal and tangent

    const global D2_short = spdiagm_const([1.0, -2.0, 1.0], [-1, 0, 1], N)/Δσ^2
    const global D2        = blkdiag(D2_short, D2_short)
    const global D2_perp  = ([[spzeros(N,N) -D2_short];[D2_short spzeros(N,N)]])

    const global D3_short = spdiagm_const([-0.5, 1.0, -1.0, 0.5], [-2, -1, 1, 2], N)/Δσ^3
    const global D3        = blkdiag(D3_short, D3_short)
    const global D3_perp  = ([[spzeros(N,N) -D3_short];[D3_short spzeros(N,N)]])

    const global D4_short = spdiagm_const([1.0, -4.0, 6.0, -4.0, 1.0], [-2, -1, 0, 1, 2], N)/Δσ^4
    const global D4        = blkdiag(D4_short, D4_short)
    const global D4_perp  = ([[spzeros(N,N) -D4_short];[D4_short spzeros(N,N)]])

    const global D5_short = spdiagm_const([-0.5, 2.0, -2.5, 2.5, -2.0, 0.5], [-3, -2, -1, 1, 2, 3], N)/Δσ^5
    const global D5        = blkdiag(D5_short, D5_short)
    const global D5_perp  = ([[spzeros(N,N) -D5_short];[D5_short spzeros(N,N)]])
end

function compute_pressure_force(coords::PointCoords, P::Params,
                                dst_f::Vector{Float64},
                                add::Bool=false)
    if add
        copy!(dst_f, dst_f - P.P*D1c_perp*vec(coords.x))
    else
        copy!(dst_f, -P.P*D1c_perp*vec(coords.x))
    end
end
function compute_pressure_force(coords::PointCoords, P::Params,
                                dst_Df::SparseMatrixCSC{Float64},
                                add::Bool=false)
    if add
        dst_Df[:] = dst_Df - P.P*D1c_perp
    else
        dst_Df[:] = -P.P*D1c_perp
    end
end

function compute_elastic_force(coords::PointCoords, coords_s::PointCoordsShifted,
                               P::Params, diffs::Differentials,
                               dst_f::Vector{Float64},
                               add::Bool=false)
    if add
        copy!(dst_f, dst_f + vec(P.K * (coords.ell.*coords.τ - coords_s.ell_m.*coords_s.τ_m)/P.Δσ))
    else
        copy!(dst_f, P.K * (coords.ell.*coords.τ - coords_s.ell_m.*coords_s.τ_m)/P.Δσ)
    end
end
function compute_elastic_force(coords::PointCoords, coords_s::PointCoordsShifted,
                               P::Params, diffs::Differentials,
                               dst_Df::SparseMatrixCSC{Float64},
                               add::Bool=false)
    if add
        dst_Df[:] = dst_Df + P.K * (AdhCommon.pointwise_projection(coords.τ)*D1p - AdhCommon.pointwise_projection(coords_s.τ_m)*D1m
                         + (AdhCommon.@bc_scalar(coords.ell).*diffs.Dτ - AdhCommon.@bc_scalar(coords_s.ell_m).*diffs.Dτ_m)
                        )/P.Δσ
    else
        dst_Df[:] = P.K * (AdhCommon.pointwise_projection(coords.τ)*D1p - AdhCommon.pointwise_projection(coords_s.τ_m)*D1m
                           + (AdhCommon.@bc_scalar(coords.ell).*diffs.Dτ - AdhCommon.@bc_scalar(coords_s.ell_m).*diffs.Dτ_m)
                          )/P.Δσ
    end
end

function compute_confinement_force(coords::PointCoords,
                                   P::Params,
                                   dst_f::Vector{Float64},
                                   add::Bool=false)
    field, ∇field, H_field = Wall.compute_field(coords.x, P; gradient=true, hessian=false)
    if add
        copy!(dst_f, dst_f + vec(-∇field))
    else
        copy!(dst,  vec(-∇field))
    end
end

function compute_confinement_force(coords::PointCoords,
                                   P::Params,
                                   dst_Df::SparseMatrixCSC{Float64},
                                   add::Bool=false)
    field, ∇field, H_field = Wall.compute_field(coords.x, P; gradient=false, hessian=true)
    if add
        dst_Df[:] = dst_Df - H_field
    else
        dst_Df[:] = -H_field
    end
end

function compute_density_increment(coords::PointCoords, coords_s::PointCoordsShifted,
                                   P::Params, F::Flags)
    """
    de/polymerization is located around the back/frontmost points, with a weight of (a, 1-2a, a)
    the total flux J = ∫(f)+ dl is fixed by the polymerization speed
    """
    a = 0.25

    x_min, x_min_idx = findmin(coords.x[:,2])
    x_max, x_max_idx = findmax(coords.x[:,2])

    # computation of the density increment
    ΔLc = 0.5*(coords.ΔL + coords_s.ΔL_m)

    front_norm = (1-2*a)*ΔLc[x_max_idx] + a*(ΔLc[idxmod(x_max_idx-1, P.N)] + ΔLc[idxmod(x_max_idx+1, P.N)])
    back_norm = (1-2*a)*ΔLc[x_min_idx] + a*(ΔLc[idxmod(x_min_idx-1, P.N)] + ΔLc[idxmod(x_min_idx+1, P.N)])

    f = zeros(P.N)

    # front section
    f[x_max_idx] = (1-2*a)*P.c/front_norm
    f[idxmod(x_max_idx-1, P.N)] = a*P.c/front_norm
    f[idxmod(x_max_idx+1, P.N)] = a*P.c/front_norm

    # back section
    f[x_min_idx] = -(1-2*a)*P.c/back_norm
    f[idxmod(x_min_idx-1, P.N)] = -a*P.c/back_norm
    f[idxmod(x_min_idx+1, P.N)] = -a*P.c/back_norm

    return f
end

function integrate(f::Array{Float64}, P::Params)
    return P.Δσ*sum(f, 1)
end

function cumsum_zero(f::Array{Float64})
    z = zeros(f)
    cumsum!(z, f, 1)
    z = 0.5z + 0.5circshift(z, 1)
    return z
end

function split_cumsum(a::Array{Float64,1}, idx)
    b = zeros(a)
    N = length(a)
    b[idx+1:end]  = cumsum(a[idx:end-1])
    b[idx-1:-1:1] = cumsum(a[idx-1:-1:1])

    m = 1
    p = N

    while abs(b[m] - b[p]) > a[p]
        if b[m] < b[p]
            b[p] = b[m] + a[p]
            m = idxmod(m+1, N)
            p = idxmod(p+1, N)
        else
            b[m] = b[p] + a[p]
            m = idxmod(m-1, N)
            p = idxmod(p-1, N)
        end
    end

    return b
end

function compute_drag_mask(coords::PointCoords, half_w::Float64, pow::Float64=2.0)
    x_min, x_min_idx = findmin(coords.x[:,2])
    x_max, x_max_idx = findmax(coords.x[:,2])

    dst_from_min = split_cumsum(coords.ΔL, x_min_idx)
    dst_from_max = split_cumsum(coords.ΔL, x_max_idx)

    profile_min = exp.(-(dst_from_min/half_w).^pow)
    profile_min .*= profile_min .> 1e-5
    profile_min .*= 0.5/sum(profile_min)
    profile_max = exp.(-(dst_from_max/half_w).^pow)
    profile_max .*= profile_max .> 1e-5
    profile_max .*= 0.5/sum(profile_max)

    return profile_min + profile_max
end

function compute_transport_force(coords::PointCoords, coords_s::PointCoordsShifted,
                                 P::Params, F::Flags,
                                 dst_f::Vector{Float64},
                                 add::Bool=false)

    f = compute_density_increment(coords, coords_s, P, F)

    # computation of the transport term
    cum_added_mass = 0.5 * cumsum_zero(f.*(coords.ΔL + coords_s.ΔL_m))

    transport_force = -vec((-0.5*P.c + cum_added_mass) .* coords.Δ2x/2P.Δσ)

    drag_mask = compute_drag_mask(coords, 1.0)
    drag_force = -vec(0.5drag_mask.*sum(f.*(coords.ΔL + coords_s.ΔL_m).*coords.x, 1))/P.Δσ

    if add
        copy!(dst_f, dst_f + transport_force + drag_force)
    else
        copy!(dst_f, transport_force + drag_force)
    end
end

function compute_transport_force(coords::PointCoords, coords_s::PointCoordsShifted,
                                 P::Params, F::Flags,
                                 diffs::Differentials,
                                 dst_Df::SparseMatrixCSC{Float64},
                                 add::Bool=false)
    f = compute_density_increment(coords, coords_s, P, F)

    # computation of the transport term
    cum_added_mass = 0.5 * cumsum_zero(f.*(coords.ΔL + coords_s.ΔL_m))

    D_transport_force = AdhCommon.@bc_scalar(-(-0.5P.c + cum_added_mass)) .* D1c

    if add
        dst_Df[:] = dst_Df + D_transport_force
    else
        dst_Df[:] = D_transport_force
    end
end

function compute_viscosity_force(coords::PointCoords, coords_s::PointCoordsShifted,
                                 P::Params, F::Flags,
                                 dst::Vector{Float64}, dst_Df::SparseMatrixCSC{Float64})
    throw("Not implemented")
    x_min, x_min_idx = findmin(view(coords.x,:,2))
    x_max, x_max_idx = findmax(view(coords.x,:,2))

    trsp_dir = 0.5P.c .* ((x_max_idx .< 1:params.N .< x_min_idx) - ((1:params.N .< x_max_idx) + (x_min_idx .< 1:params.N)))

    dst[:] = P.Ka * (1/params.δt * ((vd2τ - vd2τ_prev) .* dvΔL + (d1vd2τ - d1vd2τ_prev) .* (2params.Δσ*coords.vc)./coords.Δ2L)
                     + !F.initializing * trsp_dir .* (2P.Δσ*d2vd2τ.*coords.vc./coords.Δ2L + d1vd2τ.*dvΔL))
    # dst_Df[:] =
end

function compute_residuals(x::Vector{Float64},
                           coords::PointCoords, coords_s::PointCoordsShifted,
                           inner_coords::PointCoords, inner_coords_s::PointCoordsShifted,
                           P::Params, F::Flags,
                           dst::Vector{Float64})

    update_coords(inner_coords, P, reshape(x, (P.N, 2)))
    differentials = new_Differentials(inner_coords, inner_coords_s)

    fill!(dst, 0.0)

    compute_pressure_force(inner_coords, P, dst, true)
    compute_elastic_force(inner_coords, inner_coords_s, P, differentials, dst, true)
    if F.confine
        compute_confinement_force(inner_coords, P, dst, true)
    end
    # compute_drag_force(inner_coords, inner_coords_s, P, F, dst, true)
    compute_transport_force(inner_coords, inner_coords_s, P, F, dst, true)

    if F.innerloop
        dst[:] = vec(x) - vec(coords.x) - P.δt*dst
    else
        dst[:] = -P.δt*dst
    end

end

function compute_residuals_J(x::Vector{Float64},
                             coords::PointCoords, coords_s::PointCoordsShifted,
                             inner_coords::PointCoords, inner_coords_s::PointCoordsShifted,
                             P::Params, F::Flags,
                             dst_Df::SparseMatrixCSC{Float64})

    update_coords(inner_coords, P, reshape(x, (P.N, 2)))
    differentials = new_Differentials(inner_coords, inner_coords_s)

    if true # reset jacobian
        fill!(dst_Df.colptr, 1)
        empty!(dst_Df.rowval)
        empty!(dst_Df.nzval)
    end
    compute_pressure_force(inner_coords, P, dst_Df, true)
    compute_elastic_force(inner_coords, inner_coords_s, P, differentials, dst_Df, true)
    if F.confine
        compute_confinement_force(inner_coords, P, dst_Df, true)
    end
    compute_transport_force(inner_coords, inner_coords_s, P, F, differentials, dst_Df, true)

    dst_Df[:] = speye(2P.N) - P.δt*dst_Df
end

function wrap_residuals(coords::PointCoords, coords_s::PointCoordsShifted,
                        P::Params, F::Flags)
    inner_coords = deepcopy(coords)
    inner_coords_s = new_PointCoordsShifted(inner_coords)
    resi(x::Vector{Float64}, dst::Vector{Float64}) = compute_residuals(x, coords, coords_s, inner_coords, inner_coords_s, P, F, dst)
    resi_J(x::Vector{Float64}, dst_Df::SparseMatrixCSC{Float64}) = compute_residuals_J(x, coords, coords_s, inner_coords, inner_coords_s, P, F, dst_Df)

    return resi, resi_J
end

end # module
