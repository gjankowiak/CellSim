module Forces

using CellSimCommon
import Wall, Masks
import Utils: spdiagm_const, idxmod

struct PointCoords
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

    # centrosome
    centro_x::Vector{Float64}
    centro_angle::Vector{Float64}
end


struct PointCoordsShifted
    Δx_m::SubArray
    Δx_p::SubArray
    ΔL_m::SubArray
    ΔL_p::SubArray
    ell_m::SubArray
    τ_m::SubArray
    v_m::SubArray
    v_p::SubArray
end

mutable struct Differentials
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
         zeros(N), # d2vd2τ

         zeros(2), # centro_x
         zeros(1)) # centro_angle

    coords_shifted = new_PointCoordsShifted(coords)
    return coords, coords_shifted
end

function new_PointCoords(x::Matrix{Float64}, P::Params)
    N = size(x, 1)
    coords, coords_shifted = new_PointCoords(N)
    update_coords(coords, P, x)
    coords.centro_x[:] = sum(x, 1)/size(x,1)
    return coords, coords_shifted
end

function new_PointCoordsShifted(coords::PointCoords)
    N = size(coords.x, 1)

    return PointCoordsShifted(
        view(coords.Δx, CellSimCommon.circ_idx_m1, :),
        view(coords.Δx, CellSimCommon.circ_idx_p1, :),
        view(coords.ΔL, CellSimCommon.circ_idx_m1),
        view(coords.ΔL, CellSimCommon.circ_idx_p1),
        view(coords.ell, CellSimCommon.circ_idx_m1),
        view(coords.τ, CellSimCommon.circ_idx_m1, :),
        view(coords.v, CellSimCommon.circ_idx_m1, :),
        view(coords.v, CellSimCommon.circ_idx_p1, :))
end

function new_Differentials(coords::PointCoords, coords_s::PointCoordsShifted)
    Dτ = CellSimCommon.@bc_scalar(1./coords.ΔL).*(CellSimCommon.pointwise_projection(coords.v))*D1p_unorm
    Dτc = CellSimCommon.@bc_scalar(1./coords.Δ2L).*CellSimCommon.pointwise_projection(coords.vc)*D1c_unorm
    Dv = -CellSimCommon.@bc_scalar(1./coords.ΔL).*(CellSimCommon.pointwise_projection(coords.τ))*D1p_perp_unorm
    Dvc = -CellSimCommon.@bc_scalar(1./coords.Δ2L).*CellSimCommon.pointwise_projection(coords.τc)*D1c_perp_unorm
    return Differentials(
        Dτ, Dτc, Dv, Dvc, M_cs_minus*Dτ
    )
end

function update_coords(coords::PointCoords, P::Params, x::Matrix{Float64})
    coords.x[:] = x
    CellSimCommon.@delta!(x, coords.Δx)
    coords.Δ2x[:] = x[CellSimCommon.circ_idx_p1,:] - x[CellSimCommon.circ_idx_m1,:]
    coords.Δ2x_perp[:,1] = -coords.Δ2x[:,2]
    coords.Δ2x_perp[:,2] =  coords.Δ2x[:,1]

    coords.ΔL[:] = CellSimCommon.@entry_norm(coords.Δx)
    coords.Δ2L[:] = CellSimCommon.@entry_norm(coords.Δ2x)

    coords.ell[:] = coords.ΔL/P.Δσ - 1

    # unit tangent
    coords.τ[:] = coords.Δx./coords.ΔL

    # centered unit tangent
    coords.τc[:] = coords.Δ2x./CellSimCommon.@entry_norm(coords.Δ2x)

    # unit normal
    coords.v[:] = hcat(view(coords.τ,:,2), -view(coords.τ,:,1))

    # centered unit tangent
    coords.vc[:] = -coords.Δ2x_perp./CellSimCommon.@entry_norm(coords.Δ2x)

    # derivatives of tangent and normal unit vectors
    coords.dτc[:] = D1c*vec(coords.τc)
    coords.d2τ[:] = D2*vec(coords.τc)
    coords.d3τ[:] = D3*vec(coords.τc)

    coords.dvc[:] = D1c*vec(coords.vc)

    # auxiliary quantities
    coords.dvΔL[:] = 2P.Δσ*D1c_short*(coords.vc./coords.Δ2L)
    coords.vd2τ[:] = CellSimCommon.@dotprod(coords.vc, coords.d2τ)
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

function compute_front_back(coords::PointCoords, P::Params, F::Flags)
    if !F.circular_wall
        x_min, x_min_idx = findmin(coords.x[:,2])
        x_max, x_max_idx = findmax(coords.x[:,2])
        return ((x_min, x_min_idx), (x_max, x_max_idx))
    else
        ang = angle.(complex(coords.x[:,1] + coords.x[:,2]im))
        x_min, x_min_idx = findmin(ang)
        x_max, x_max_idx = findmax(ang)
        if (x_max - x_min) > pi
            x_min, x_max = x_max, x_min
            x_min_idx, x_max_idx = x_max_idx, x_min_idx
        end
        return ((x_min, x_min_idx), (x_max, x_max_idx))
    end
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
        dst_Df[:] = dst_Df + P.K * (CellSimCommon.pointwise_projection(coords.τ)*D1p - CellSimCommon.pointwise_projection(coords_s.τ_m)*D1m
                         + (CellSimCommon.@bc_scalar(coords.ell).*diffs.Dτ - CellSimCommon.@bc_scalar(coords_s.ell_m).*diffs.Dτ_m)
                        )/P.Δσ
    else
        dst_Df[:] = P.K * (CellSimCommon.pointwise_projection(coords.τ)*D1p - CellSimCommon.pointwise_projection(coords_s.τ_m)*D1m
                           + (CellSimCommon.@bc_scalar(coords.ell).*diffs.Dτ - CellSimCommon.@bc_scalar(coords_s.ell_m).*diffs.Dτ_m)
                          )/P.Δσ
    end
end

function compute_confinement_force(coords::PointCoords,
                                   P::Params, F::Flags, plotables::CellSimCommon.Plotables,
                                   dst_f::Vector{Float64},
                                   add::Bool=false, weighted::Bool=false)
    field, ∇field, H_field = Wall.compute_field(coords.x, P, F; gradient=true, hessian=false)
    if weighted
        copy!(∇field, coords.Δ2L.*∇field/2P.Δσ)
    end
    if add
        copy!(dst_f, dst_f + vec(-∇field))
    else
        copy!(dst,  vec(-∇field))
    end
    plotables.field[:] = field
    plotables.∇field[:] = ∇field
end

function compute_confinement_force(coords::PointCoords,
                                   P::Params, F::Flags,
                                   dst_Df::SparseMatrixCSC{Float64},
                                   add::Bool=false, weighted::Bool=false)
    field, ∇field, H_field = Wall.compute_field(coords.x, P, F; gradient=true, hessian=true)
    if weighted
        aux = ([[D1c_short_unorm.*coords.τc[:,1] D1c_short_unorm.*coords.τc[:,2]]
                [D1c_short_unorm.*coords.τc[:,1] D1c_short_unorm.*coords.τc[:,2]]])
        H_field = (CellSimCommon.@bc_scalar(coords.Δ2L).*H_field + vec(∇field).*aux)/2P.Δσ
    end
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

    ((x_min, x_min_idx), (x_max, x_max_idx)) = compute_front_back(coords, P, F)

    # computation of the density increment
    ΔLc = 0.5*(coords.ΔL + coords_s.ΔL_m)

    mask_front, mask_back = compute_mask(coords, P, F, P.mass_gauss_width, P.mass_gauss_power)

    front_norm = sum(mask_front.*ΔLc)
    back_norm = sum(mask_back.*ΔLc)

    f = zeros(P.N)
    f = P.c*(mask_front/front_norm - mask_back/back_norm)

    return f
end

function integrate(f::Array{Float64}, P::Params)
    return P.Δσ*sum(f, 1)
end

function cumsum_zero(f::Array{Float64})
    z = zeros(f)
    cumsum!(z, f, 1)
    z -= f
    # z = 0.5z + 0.5circshift(z, 1)
    return z
end

function split_cumsum(a::Array{Float64,1}, idx)
    N = length(a)

    aux = zeros(a)

    dist_left = copy(a)
    dist_left[idx] = 0.0
    circshift!(aux, dist_left, N-idx)
    reverse!(aux)
    cumsum!(aux, aux)
    reverse!(aux)
    circshift!(dist_left, aux, idx-N)

    dist_right = copy(a)
    dist_right[idx] = 0.0
    circshift!(aux, dist_right, 1-idx)
    cumsum!(aux, aux)
    circshift!(dist_right, aux, idx-1)

    return min.(dist_left, dist_right)
end

function compute_mask(coords::PointCoords, P::Params, F::Flags, half_w::Float64, pow::Float64=2.0)
    ((x_min, x_min_idx), (x_max, x_max_idx)) = compute_front_back(coords, P, F)

    dst_from_min = split_cumsum(coords.ΔL, x_min_idx)
    dst_from_max = split_cumsum(coords.ΔL, x_max_idx)

    profile_min = exp.(-(dst_from_min/half_w).^pow)
    profile_min .*= profile_min .> 1e-5
    profile_min .*= 0.5/sum(profile_min)
    profile_max = exp.(-(dst_from_max/half_w).^pow)
    profile_max .*= profile_max .> 1e-5
    profile_max .*= 0.5/sum(profile_max)

    return profile_max, profile_min
end

function compute_transport_force(coords::PointCoords, coords_s::PointCoordsShifted,
                                 P::Params, F::Flags, plt::CellSimCommon.Plotables,
                                 dst_f::Vector{Float64},
                                 add::Bool=false)

    f = compute_density_increment(coords, coords_s, P, F)
    plt.mass_source[:] = 0.5(f.*(coords.ΔL + coords_s.ΔL_m))

    # computation of the transport term
    plt.mass_source_int[:] = 0.5cumsum_zero(f.*(coords.ΔL + coords_s.ΔL_m))

    plt.transport_force[:] = (P.Δσ*sum(plt.mass_source_int)-plt.mass_source_int).*coords.Δ2x/2P.Δσ

    drag_mask_f, drag_mask_b = compute_mask(coords, P, F, P.drag_gauss_width, P.drag_gauss_power)
    drag_mask = drag_mask_f + drag_mask_b

    # caution, the centered difference operator is not skewsymmetric
    plt.drag_force[:] = drag_mask.*sum(plt.mass_source_int .* coords.Δ2x/2P.Δσ, 1)
    # plt.drag_force[:] = drag_mask.*sum(plt.mass_source_int .* coords.Δx/P.Δσ, 1)

    if add
        copy!(dst_f, dst_f + vec(plt.transport_force + plt.drag_force))
    else
        copy!(dst_f, vec(plt.transport_force + plt.drag_force))
    end
end

function compute_transport_force(coords::PointCoords, coords_s::PointCoordsShifted,
                                 P::Params, F::Flags, plt::CellSimCommon.Plotables,
                                 diffs::Differentials,
                                 dst_Df::SparseMatrixCSC{Float64},
                                 add::Bool=false)

    D_transport_force = CellSimCommon.@bc_scalar(P.Δσ*sum(plt.mass_source_int)-plt.mass_source_int) .* D1c
    # D_transport_force = CellSimCommon.@bc_scalar(P.Δσ*sum(plt.mass_source_int)-plt.mass_source_int) .* D1p

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
    ((x_min, x_min_idx), (x_max, x_max_idx)) = compute_front_back(coords, P, F)

    trsp_dir = 0.5P.c .* ((x_max_idx .< 1:params.N .< x_min_idx) - ((1:params.N .< x_max_idx) + (x_min_idx .< 1:params.N)))

    dst[:] = P.Ka * (1/params.δt * ((vd2τ - vd2τ_prev) .* dvΔL + (d1vd2τ - d1vd2τ_prev) .* (2params.Δσ*coords.vc)./coords.Δ2L)
                     + !F.initializing * trsp_dir .* (2P.Δσ*d2vd2τ.*coords.vc./coords.Δ2L + d1vd2τ.*dvΔL))
    # dst_Df[:] =
end

function compute_residuals(x::Vector{Float64},
                           coords::PointCoords, coords_s::PointCoordsShifted,
                           inner_coords::PointCoords, inner_coords_s::PointCoordsShifted,
                           P::Params, F::Flags, plotables::CellSimCommon.Plotables,
                           dst::Vector{Float64})

    update_coords(inner_coords, P, reshape(x, (P.N, 2)))
    differentials = new_Differentials(inner_coords, inner_coords_s)

    fill!(dst, 0.0)

    compute_pressure_force(inner_coords, P, dst, true)
    compute_elastic_force(inner_coords, inner_coords_s, P, differentials, dst, true)
    if F.confine
        compute_confinement_force(inner_coords, P, F, plotables, dst, true, F.weighted_confinement)
    end
    if F.polymerize
        compute_transport_force(inner_coords, inner_coords_s, P, F, plotables, dst, true)
    end

    if F.innerloop
        dst[:] = vec(x) - vec(coords.x) - P.δt*dst
    else
        dst[:] = -P.δt*dst
    end

end

function compute_residuals_J(x::Vector{Float64},
                             coords::PointCoords, coords_s::PointCoordsShifted,
                             inner_coords::PointCoords, inner_coords_s::PointCoordsShifted,
                             P::Params, F::Flags, plotables::CellSimCommon.Plotables,
                             dst_Df::SparseMatrixCSC{Float64})

    update_coords(inner_coords, P, reshape(x, (P.N, 2)))
    differentials = new_Differentials(inner_coords, inner_coords_s)

    if true # reset jacobian
        fill!(dst_Df.colptr, 1)
        empty!(dst_Df.rowval)
        empty!(dst_Df.nzval)
    end
    compute_pressure_force(inner_coords, P, dst_Df, true)
    compute_elastic_force(inner_coords, inner_coords_s,
                          P, differentials, dst_Df, true)
    if F.confine
        compute_confinement_force(inner_coords, P, F, dst_Df, true, F.weighted_confinement)
    end

    if F.polymerize
        compute_transport_force(inner_coords, inner_coords_s,
                                P, F, plotables,
                                differentials, dst_Df, true)
    end

    dst_Df[:] = speye(2P.N) - P.δt*dst_Df
end

function wrap_residuals(coords::PointCoords, coords_s::PointCoordsShifted,
                        P::Params, F::Flags, plotables::CellSimCommon.Plotables)
    inner_coords = deepcopy(coords)
    inner_coords_s = new_PointCoordsShifted(inner_coords)
    resi(x::Vector{Float64}, dst::Vector{Float64}) =
        compute_residuals(x, coords, coords_s, inner_coords, inner_coords_s, P, F, plotables, dst)
    resi_J(x::Vector{Float64}, dst_Df::SparseMatrixCSC{Float64}) =
        compute_residuals_J(x, coords, coords_s, inner_coords, inner_coords_s, P, F, plotables, dst_Df)

    return resi, resi_J
end

end # module
