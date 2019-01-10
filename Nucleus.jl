module Nucleus

const DEBUG = "DEBUG" in keys(ENV)

import Cortex: PointCoords, PointCoordsShifted

using CellSimCommon
const CSC = CellSimCommon

import SparseArrays
const SA = SparseArrays

import LinearAlgebra
const LA = LinearAlgebra

struct NucleusCoords
    # node coordinates
    Y::Matrix{Float64}

    # tangential velocity
    α::Vector{Float64}

    # normal velocity
    β::Vector{Float64}

    # curvature
    k::Vector{Float64}

    # finite volume size, i.e. |Y_i - Y_i-1|
    r::Vector{Float64}

    # dual finite volume size, i.e. 0.5(r_i+r_i+1)
    q::Vector{Float64}

    # tangential angle
    θ::Vector{Float64}

    # normal vector
    n::Matrix{Float64}

    # local length element
    η::Vector{Float64}

    # nuclear membrane length
    L::Float64

    circ_idx::CSC.CircIdx
end

# TODO
# - introduce views for shifted arrays
# - FIX the handl

function g(x::Vector{Float64}, α::Float64)
    return -min.(α*x.-1, 0.0).^2 .* log.(α*x)
end

function g_p(x::Vector{Float64}, α::Float64)
    return -(2α*min.(α*x.-1, 0.0).*log.(α*x) .+ min.(α*x.-1, 0.0).^2 ./x)
end

function g_pp(x::Vector{Float64}, α::Float64)
    return -(2α^2*(x.<(1/α)).*log.(α*x) .+ 4α*min.(α*x.-1, 0.0)./x .- min.(α*x.-1, 0.0).^2 ./x.^2)
end

function compute_contact_force(pots::CSC.InteractionPotentials,
                               cor_coords::PointCoords, c::NucleusCoords,
                               P::Params, F::Flags)
    for i in 1:P.Nnuc
        xy = cor_coords.x .- c.Y[i,:]'
        d = sqrt.(sum(abs2, xy; dims=2))
        xy_norm = xy ./ d

        pot = g(vec(d), 10.0)
        ∇pot = -xy_norm .* g_p(vec(d), 10.0)

        pots.C_∇W[:] = pots.C_∇W - ∇pot
        pots.N_W[i] = pots.N_W[i] + sum(pot)
        pots.N_∇W[i,:] = pots.N_∇W[i,:] + vec(sum(∇pot; dims=1))
    end
end

function compute_centronuclear_force(pots::CSC.InteractionPotentials,
                                     cor_coords::PointCoords, c::NucleusCoords,
                                     P::Params, F::Flags)

    Yc = c.Y .- cor_coords.centro_x'
    d = sqrt.(sum(abs2, Yc; dims=2))
    Yc_norm = Yc ./ d

    pots.CS_∇W[:] = pots.CS_∇W + sum(P.N_kc * (d .- P.N_l0c) .* Yc_norm; dims=1)

    pots.N_W[:] = pots.N_W + 0.5*P.N_kc/P.Nnuc*(d .- P.N_l0c).^2
    pots.N_∇W[:] = pots.N_∇W - repeat(pots.CS_∇W, P.Nnuc)/P.Nnuc
end

function update_alphabeta(c::NucleusCoords, new_c::NucleusCoords,
                          W::Vector{Float64}, ∇W::Matrix{Float64},
                          P::Params, F::Flags)
    circ_idx = new_c.circ_idx

    c.β[:] = (
              P.N_kb./c.r .*((c.k[circ_idx.p1]-c.k)./c.q - (c.k - c.k[circ_idx.m1])./c.q[circ_idx.m1])
              + 0.5P.N_kb*c.k.^3
              .- P.N_P
              + 0.5*(W+W[circ_idx.m1]).*c.k
              - 0.5vec(sum(c.n.*(∇W+∇W[circ_idx.m1,:]); dims=2))
             )

    B = sum(c.r .* c.k .* c.β) / c.L
    new_c.α[1] = 0
    for i in 2:P.Nnuc
        new_c.α[i] = (new_c.α[i-1]
                           + c.r[i]*(-c.k[i]*c.β[i]
                                     +B)
                           + (c.L/P.Nnuc - c.r[i])*P.N_ω)
    end
end

function update_r(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags)
    new_c.η[:] = c.η + P.δt*(c.k .* c.β
                                               +(new_c.α - new_c.α[new_c.circ_idx.m1])
                                               ./c.r)
    new_c.r[:] = exp.(new_c.η)
end

function finite_differences_3(new_c::NucleusCoords,
                              Dm2::Vector{Float64}, Dm1::Vector{Float64},
                              D0::Vector{Float64},
                              Dp1::Vector{Float64}, Dp2::Vector{Float64};
                              prefactor::Float64=1.0)
    c = new_c
    circ_idx = c.circ_idx
    Dm2[:] = @. prefactor/ (c.q[circ_idx.m1] * c.q[circ_idx.m2] * c.r[circ_idx.m1])
    Dp2[:] = @. prefactor/ (c.q * c.q[circ_idx.p1] * c.r[circ_idx.p1])

    Dm1[:] = @. -prefactor * (
             1 / (c.r * c.q * c.q[circ_idx.m1])
            +1 / (c.r * c.q[circ_idx.m1]^2)
            +1 / (c.r[circ_idx.m1] * c.q[circ_idx.m1]^2)
            +1 / (c.r[circ_idx.m1] * c.q[circ_idx.m1] * c.q[circ_idx.m2])
           )
    Dm1[:] = @. -prefactor * (
             1 / (c.r[circ_idx.p1] * c.q[circ_idx.p1] * c.q)
            +1 / (c.r[circ_idx.p1] * c.q^2)
            +1 / (c.r * c.q^2)
            +1 / (c.r * c.q * c.q[circ_idx.m1])
           )
    D0[:] = -(Dm2 + Dm1 + Dp1 + Dp2)
end

function spdiagm_wrap(N::Int64, kv::Pair{Int64, Vector{Float64}}...)
    new_kv = ()
    for (k, v) in kv
        @assert length(v) == N
        if k == 0
            new_kv = (new_kv..., 0 => v)
        else
            new_kv = (new_kv..., k => v[1:N-abs(k)])
            new_kv = (new_kv..., -sign(k)*(N-abs(k)) => v[N-abs(k)+1:end])
        end
    end
    return SA.spdiagm(new_kv...)
end

function update_K(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::CSC.TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, f) = CSC.@ta6_tuple(temparrays)

    finite_differences_3(new_c, Dm2, Dm2, D0, Dp1, Dp2; prefactor=P.N_kb)

    Dm1[:] = Dm1 + @. 0.5 * (new_c.α[circ_idx.m1] - (W[circ_idx.m1] + W[circ_idx.m2])./new_c.q[circ_idx.m1])
    Dp1[:] = Dp1 + @. 0.5 * (-new_c.α - (W[circ_idx.p1] + W)./new_c.q)

    D0[:] = D0 +(
           new_c.r./P.δt
           + 0.5*(new_c.α - new_c.α[circ_idx.m1])
           + 0.5*(W + W[circ_idx.m1]).*(1 ./new_c.q + 1 ./new_c.q[circ_idx.m1])
           + c.r .* c.k .* c.β)

    f[:] = (
            new_c.r.*c.k/P.δt
            - 0.5P.N_kb*(c.k[circ_idx.p1].^3 - c.k.^3)./c.q
            + 0.5P.N_kb*(c.k.^3 - c.k[circ_idx.m1].^3)./c.q[circ_idx.p1]
            + 0.5*vec( CSC.@dotprod(∇W[circ_idx.p1,:] + ∇W, c.n[circ_idx.p1,:]) ./ c.q
                   +CSC.@dotprod(∇W + ∇W[circ_idx.m1,:], c.n) .* (1 ./c.q + 1 ./c.q[circ_idx.m1])
                   +CSC.@dotprod(∇W[circ_idx.m1,:]+∇W[circ_idx.m2,:], c.n[circ_idx.m1,:]) ./ c.q[circ_idx.m1]
                  )
           )

    M = spdiagm_wrap(P.Nnuc,
                     -2 => Dm2,
                     -1 => Dm1,
                     0 => D0,
                     1 => Dp1,
                     2 => Dp2)

    display(Matrix(M))
    new_c.k[:] = M\f
end

function update_θ(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::CSC.TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, f) = CSC.@ta6_tuple(temparrays)

    finite_differences_3(new_c, Dm2, Dm2, D0, Dp1, Dp2; prefactor=P.N_kb)

    Dm1[:] = Dm1 + 0.5new_c.α[circ_idx.m1] - W[circ_idx.m1]./new_c.q[circ_idx.m1]
    Dp1[:] = Dp1 - 0.5new_c.α - W./new_c.q

    D0[:] = D0 + (
                  new_c.r/P.δt
                  + 0.5*(new_c.α + new_c.α[circ_idx.m1])
                  - (W./new_c.q + W[circ_idx.m1]./new_c.q[circ_idx.m1])
                 )

    f[:] = (new_c.r.*c.θ/P.δt
            - 0.5P.N_kb*(((c.θ[circ_idx.p1] - c.θ)./c.q).^3 + ((c.θ - c.θ[circ_idx.m1])./c.q[circ_idx.m1]).^3)
            + 0.5*vec(CSC.@dotprod(∇W, (c.n[circ_idx.p1,:] + c.n)) + CSC.@dotprod(∇W[circ_idx.m1,:], (c.n + c.n[circ_idx.m1,:])))
           )

    M = spdiagm_wrap(P.Nnuc,
                     -2 => Dm2,
                     -1 => Dm1,
                     0 => D0,
                     1 => Dp1,
                     2 => Dp2)

    new_c.θ[:] = M\f

    new_c.n[:] = [sin.(new_c.θ) -cos.(new_c.θ)]
end

function update_Y(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::CSC.TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, _) = CSC.@ta6_tuple(temparrays)

    finite_differences_3(new_c, Dm2, Dm2, D0, Dp1, Dp2; prefactor=P.N_kb)

    Dm1[:] = Dm1 + 1.5P.N_kb*new_c.k.^2 ./new_c.r + 0.5new_c.α - W./new_c.r
    Dp1[:] = Dp1 + 1.5P.N_kb*new_c.k[circ_idx.p1].^2 ./new_c.r[circ_idx.p1] - 0.5new_c.α - W./new_c.r[circ_idx.p1]

    D0[:] = new_c.q/P.δt -(Dm2 + Dm1 + Dp1 + Dp2)

    n = new_c.n
    n_p1 = circshift(n, -1)

    f = (
         vec(new_c.q .* c.Y)/P.δt
         -0.5P.N_P*vec(new_c.r[circ_idx.p1].*n_p1 + new_c.r.*n)
         -0.5*vec(
                c.r[circ_idx.p1] .* vec(CSC.@dotprod(c.n[circ_idx.p1,:], ∇W)) .* c.n[circ_idx.p1,:]
               +c.r .* vec(CSC.@dotprod(c.n, ∇W)) .* c.n
              )
        )


    _M = spdiagm_wrap(P.Nnuc,
                      -2 => Dm2,
                      -1 => Dm1,
                      0 => D0,
                      1 => Dp1,
                      2 => Dp2)

    M = CSC.@repdiagblk(_M, 2)

    new_c.Y[:] = M\f
end

function initialize_coords(P::Params, F::Flags, cortex_c::PointCoords)
    # The nucleus is initialized as a circle centered on the
    # barycenter of the cell, with radius P.N_r_init

    bc = sum(cortex_c.x; dims=1) / P.N
    t = collect(range(0, stop=2pi, length=P.Nnuc+1))[1:P.Nnuc]

    nc = NucleusCoords(
        P.N_r_init*[cos.(t) sin.(t)] .+ bc, # Y
        zeros(P.Nnuc), # α
        zeros(P.Nnuc), # β
        1/P.N_r_init * ones(P.Nnuc), # k
        sin(π/P.Nnuc)*P.N_r_init * ones(P.Nnuc), # r
        sin(π/P.Nnuc)*P.N_r_init * ones(P.Nnuc), # q
        collect(2π/P.Nnuc * (-1:(P.Nnuc-2))), # θ
        zeros(P.Nnuc, 2), # n
        zeros(P.Nnuc), # η
        2π*P.N_r_init, # L
        CSC.init_circ_idx(P.Nnuc) # circ_idx
       )

    nc.Y[:] = nc.Y .+ bc
    nc.n[:] = [sin.(nc.θ) cos.(nc.θ)]

    return nc
end

function update_coords(c::NucleusCoords, new_c::NucleusCoords,
                       potentials::CSC.InteractionPotentials,
                       P::Params, F::Flags, temparrays::CSC.TempArrays6)

    N_W = potentials.N_W
    N_∇W = potentials.N_∇W

    update_alphabeta(c, new_c, N_W, N_∇W, P, F)
    update_r(c, new_c, N_W, N_∇W, P, F)
    update_K(c, new_c, N_W, N_∇W, P, F, temparrays)
    update_θ(c, new_c, N_W, N_∇W, P, F, temparrays)
    update_Y(c, new_c, N_W, N_∇W, P, F, temparrays)
end

end # module Nucleus
