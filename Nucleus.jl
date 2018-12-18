module Nucleus

const DEBUG = "DEBUG" in keys(ENV)

import Cortex: PointCoords, PointCoordsShifted

using CellSimCommon

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

    # nuclear membrane lenght
    L::Float64

    circ_idx::CircIdx
end

# TODO
# - introduce views for shifted arrays
# - FIX the handling of W / ∇W

function compute_contact_force(cor_coords::PointCoords, c::NucleusCoords,
                               P::Params, F::Flags)
end

function compute_centronuclear_force(cor_coords::PointCoords, c::NucleusCoords,
                                     P::Params, F::Flags)
end

function update_alphabeta(c::NucleusCoords, new_c::NucleusCoords,
                          W::Vector{Float64}, ∇W::Matrix{Float64},
                          P::Params, F::Flags)
    c.β[:] = (
              P.N_kb./c.r .*((c.k[circ_idx.p1]-c.k)./c.q - (c.k - c.k[circ_idx.m1])./c.q[circ_idx.m1])
              + 0.5P.N_kb*c.k.^3
              - P.N_P
              + W.*c.k
              - vec(sum(c.n.*∇W; dims=2))
             )

    B = sum(c.r .* c.k .* c.β) / c.L
    new_c.α[1] = 0
    for i in 2:P.Nnuc
        new_c.α[i] = (new_c.α[i-1]
                           + c.r[i]*(-c.k[i]*c.β[i]
                                              +c.B)
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
            new_kv = (new_kv..., k => v)
        else
            new_kv = (new_kv..., k => v[1:N-abs(k)])
            new_kv = (new_kv..., -sign(k)*(N-k) => v[N-abs(k)+1:end])
        end
    end
    println(new_kv)
    return SA.spdiagm(new_kv...)
end

function update_K(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, f) = @ta6_tuple(temparrays)

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
            + (W[circ_idx.p1] - W).*(tan.(c.θ[circ_idx.p1]) - cot.(c.θ[circ_idx.p1]))./(c.q.*c.r[circ_idx.p1])
            - (W - W[circ_idx.m1]).*(tan.(c.θ) - cot.(c.θ)).*(1 ./c.q + 1./c.q[circ_idx.m1])
            + (W[circ_idx.m1] - W[circ_idx.m2]).*(tan.(c.θ[circ_idx.m1]) - cot.(c.θ[circ_idx.m1]))./(c.q[circ_idx.m1].*c.r[circ_idx.m1])
           )

    M = spdiagm_wrap(-2 => Dm2,
                   -1 => Dm1,
                   0 => D0,
                   1 => Dp1,
                   2 => Dp2)

    new_c.k[:] = M\f
end

function update_θ(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, f) = @ta6_tuple(temparrays)

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
            + 0.5*(W[circ_idx.p1] - W).*(tan.(c.θ[circ_idx.p1]) - cot.(c.θ[circ_idx.p1]))./c.r[circ_idx.p1]
            - 0.5*(W[circ_idx.m1] - W[circ_idx.m2]).*(tan.(c.θ[circ_idx.m1]) - cot.(c.θ[circ_idx.m1]))./c.r[circ_idx.m1]
           )

    M = spdiagm_wrap(-2 => Dm2,
                   -1 => Dm1,
                   0 => D0,
                   1 => Dp1,
                   2 => Dp2)

    new_c.θ[:] = M\f

    new_c.n[:] = [sin.(new_c.θ) -cos.(new_c.θ)]
end

function update_Y(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::TempArray6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, _) = @ta6_tuple(temparrays)

    finite_differences_3(new_c, Dm2, Dm2, D0, Dp1, Dp2; prefactor=P.N_kb)

    Dm1[:] = Dm1 + 1.5P.N_kb*new_c.k.^2./new_c.r + 0.5new_c.α - W./new_c.r
    Dp1[:] = Dp1 + 1.5P.N_kb*new_c.k[circ_idx.p1].^2./new_c.r[circ_idx.p1] - 0.5new_c.α - W./new_c.r[circ_idx.p1]

    D0[:] = new_c.q/P.δt -(Dm2 + Dm1 + Dp1 + Dp2)

    n = new_c.n
    n_p1 = circshift(n, -1)

    f = (
         vec(new_c.q .* c.Y)/P.δt
         -0.5P.N_P*vec(new_c.r[circ_idx.p1].*n_p1 + new_c.r.*n)
         -0.5*vec((W[circ_idx.p1] - W).*(tan.(new_c.θ[circ_idx.p1]) - cot.(new_c.θ[circ_idx.p1])).*n_p1)
         -0.5*vec((W - W[circ_idx.m1]).*(tan.(new_c.θ) - cot.(new_c.θ)).*n)
        )

    _M = spdiagm_wrap(-2 => Dm2,
                   -1 => Dm1,
                   0 => D0,
                   1 => Dp1,
                   2 => Dp2)

    M = @repdiagblk(_M, 2)

    new_c.Y[:] = M\f
end

function update_coords(c::NucleusCoords, new_c::NucleusCoords,
                       W::Vector{Float64}, ∇W::Matrix{Float64},
                       P::Params, F::Flags)

    update_alphabeta(c, new_c, W, ∇W, P, F, temparrays)
    update_r(c, new_c, W, ∇W, P, F, temparrays)
    update_K(c, new_c, W, ∇W, P, F, temparrays)
    update_θ(c, new_c, W, ∇W, P, F, temparrays)
    update_Y(c, new_c, W, ∇W, P, F, temparrays)
end

end # module Nucleus
