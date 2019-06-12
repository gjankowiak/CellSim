module Nucleus

import Cortex: PointCoords, PointCoordsShifted

using CellSimCommon
const CSC = CellSimCommon

import SparseArrays
const SA = SparseArrays

import LinearAlgebra
const LA = LinearAlgebra

mutable struct NucleusCoords
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

function copy(dst::NucleusCoords, src::NucleusCoords)
    dst.Y[:] = src.Y
    dst.α[:] = src.α
    dst.β[:] = src.β
    dst.k[:] = src.k
    dst.r[:] = src.r
    dst.q[:] = src.q
    dst.θ[:] = src.θ
    dst.n[:] = src.n
    dst.η[:] = src.η
    dst.L = src.L
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

function compute_contact_force(pots::CSC.InteractionPotentials,
                               cor_coords::PointCoords, c::NucleusCoords,
                               P::Params, F::Flags)

    # FIXME this needs to be normalized!
    if P.N_kcont > 0
        for i in 1:P.Nnuc
            xy = cor_coords.x .- c.Y[i,:]'
            d = sqrt.(sum(abs2, xy; dims=2))
            xy_norm = xy ./ d

            pot = P.N_kcont*g(vec(d), P.N_αcont)
            ∇pot = -P.N_kcont * xy_norm .* g_p(vec(d), P.N_αcont)

            pots.C_∇W[:] = pots.C_∇W - ∇pot
            pots.N_W[i] = pots.N_W[i] + sum(pot)
            pots.N_∇W[i,:] = pots.N_∇W[i,:] + vec(sum(∇pot; dims=1))
        end
    end

    pots.N_W[:] = pots.N_W .+ P.N_W0
end

function compute_centronuclear_force(pots::CSC.InteractionPotentials,
                                     cor_coords::PointCoords, c::NucleusCoords,
                                     P::Params, F::Flags)

    if P.N_kc > 0
        Yc = c.Y .- cor_coords.centro_x'
        d = sqrt.(sum(abs2, Yc; dims=2))
        Yc_norm = Yc ./ d

        pots.CS_∇W[:] = pots.CS_∇W + sum(P.N_kc * (d .- P.N_l0c) .* Yc_norm; dims=1)

        pots.N_W[:] = pots.N_W + 0.5*P.N_kc/P.Nnuc*(d .- P.N_l0c).^2
        pots.N_∇W[:] = pots.N_∇W - repeat(pots.CS_∇W, P.Nnuc)/P.Nnuc
    end
end

function update_alphabeta(c::NucleusCoords, new_c::NucleusCoords,
                          W::Vector{Float64}, ∇W::Matrix{Float64},
                          P::Params, F::Flags)
    circ_idx = new_c.circ_idx

    c.L = sum(c.r)

    area = 0.25*sum(c.q.*(CSC.@dotprod(c.n, c.Y) + CSC.@dotprod(c.n[circ_idx.p1,:], c.Y)))
    PA = -P.N_mu*(area-P.N_target_area)

    c.β[:] = (
              P.N_kb./c.r .*( (c.k[circ_idx.p1] - c.k)./c.q
                             -(c.k - c.k[circ_idx.m1])./c.q[circ_idx.m1])
              + 0.5P.N_kb*c.k.^3
              .- (P.N_P - PA)
              - 0.5*(W+W[circ_idx.m1]).*c.k
              - 0.5vec(sum(c.n.*(∇W+∇W[circ_idx.m1,:]); dims=2))
             )

    # for plotting only
    new_c.β[:] = c.β

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
    new_c.η[:] = c.η + P.δt*(c.k .* c.β +(new_c.α - new_c.α[new_c.circ_idx.m1])./c.r)
    new_c.r[:] = exp.(new_c.η)
    new_c.q[:] = 0.5*(new_c.r + new_c.r[new_c.circ_idx.p1])
end

function finite_differences_3(c::NucleusCoords,
                              Dm2::Vector{Float64}, Dm1::Vector{Float64},
                              D0::Vector{Float64},
                              Dp1::Vector{Float64}, Dp2::Vector{Float64};
                              prefactor::Float64=1.0, dual::Bool=false)
    circ_idx = c.circ_idx
    if dual
        Dm2[:] = @. prefactor/ (c.q[circ_idx.m1] * c.q[circ_idx.m2] * c.r[circ_idx.m1])
        Dp2[:] = @. prefactor/ (c.q * c.q[circ_idx.p1] * c.r[circ_idx.p1])

        Dm1[:] = @. -prefactor * (
                 1 / (c.r * c.q * c.q[circ_idx.m1])
                +1 / (c.r * c.q[circ_idx.m1]^2)
                +1 / (c.r[circ_idx.m1] * c.q[circ_idx.m1]^2)
                +1 / (c.r[circ_idx.m1] * c.q[circ_idx.m1] * c.q[circ_idx.m2])
               )
        Dp1[:] = @. -prefactor * (
                 1 / (c.r[circ_idx.p1] * c.q[circ_idx.p1] * c.q)
                +1 / (c.r[circ_idx.p1] * c.q^2)
                +1 / (c.r * c.q^2)
                +1 / (c.r * c.q * c.q[circ_idx.m1])
               )
    else
        Dm2[:] = @. prefactor/ (c.r * c.r[circ_idx.m1] * c.q[circ_idx.m1])
        Dp2[:] = @. prefactor/ (c.r[circ_idx.p1] * c.r[circ_idx.p2] * c.q[circ_idx.p1])

        Dm1[:] = @. -prefactor * (
                 1 / (c.q[circ_idx.m1] * c.r[circ_idx.m1] * c.r)
                +1 / (c.q[circ_idx.m1] * c.r^2)
                +1 / (c.q * c.r^2)
                +1 / (c.q * c.r[circ_idx.p1] * c.r)
               )
        Dp1[:] = @. -prefactor * (
                 1 / (c.q * c.r * c.r[circ_idx.p1])
                +1 / (c.q * c.r[circ_idx.p1]^2)
                +1 / (c.q[circ_idx.p1] * c.r[circ_idx.p1]^2)
                +1 / (c.q[circ_idx.p1] * c.r[circ_idx.p1] * c.r[circ_idx.p2])
               )
    end
    D0[:] = -(Dm2 + Dm1 + Dp1 + Dp2)
end

function spdiagm_wrap(N::Int64, kv::Pair{Int64, Vector{Float64}}...)
    new_kv = ()
    for (k, v) in kv
        @assert length(v) == N
        if k == 0
            new_kv = (new_kv..., 0 => v)
        else
            if k > 0
                new_kv = (new_kv..., k => v[1:N-k])
                new_kv = (new_kv..., k-N => v[N-k+1:end])
            else
                new_kv = (new_kv..., k => v[-k+1:end])
                new_kv = (new_kv..., (N+k) => v[1:-k])
            end
        end
    end
    return SA.spdiagm(new_kv...)
end

function update_K(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::CSC.TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, f) = CSC.@ta6_tuple(temparrays)

    finite_differences_3(new_c, Dm2, Dm1, D0, Dp1, Dp2; prefactor=P.N_kb)

    Dm1[:] = Dm1 + @. 0.5 * (new_c.α[circ_idx.m1] - (W[circ_idx.m1] + W[circ_idx.m2])./new_c.q[circ_idx.m1])
    Dp1[:] = Dp1 + @. 0.5 * (-new_c.α - (W[circ_idx.p1] + W)./new_c.q)

    D0[:] = D0 +(
           new_c.r./P.δt
           + 0.5*(new_c.α - new_c.α[circ_idx.m1])
           + 0.5*(W + W[circ_idx.m1]).*(1 ./new_c.q + 1 ./new_c.q[circ_idx.m1])
           + new_c.r .* c.k .* c.β)

    f[:] = (
            new_c.r.*c.k/P.δt
            - 0.5P.N_kb*(c.k[circ_idx.p1].^3 - c.k.^3)./new_c.q
            + 0.5P.N_kb*(c.k.^3 - c.k[circ_idx.m1].^3)./new_c.q[circ_idx.m1]
            + 0.5*vec( CSC.@dotprod(∇W[circ_idx.p1,:] + ∇W, c.n[circ_idx.p1,:]) ./ new_c.q
                   -CSC.@dotprod(∇W + ∇W[circ_idx.m1,:], new_c.n) .* (1 ./new_c.q + 1 ./new_c.q[circ_idx.m1])
                   +CSC.@dotprod(∇W[circ_idx.m1,:]+∇W[circ_idx.m2,:], new_c.n[circ_idx.m1,:]) ./ new_c.q[circ_idx.m1]
                  )
           )

    M = spdiagm_wrap(P.Nnuc,
                     -2 => Dm2,
                     -1 => Dm1,
                     0 => D0,
                     1 => Dp1,
                     2 => Dp2)

    new_c.k[:] = M\f
end

function update_θ(c::NucleusCoords, new_c::NucleusCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags, temparrays::CSC.TempArrays6)
    circ_idx = new_c.circ_idx

    (Dm2, Dm1, D0, Dp1, Dp2, f) = CSC.@ta6_tuple(temparrays)

    finite_differences_3(new_c, Dm2, Dm1, D0, Dp1, Dp2; prefactor=P.N_kb)

    Dm1[:] = Dm1 + 0.5new_c.α[circ_idx.m1] - W[circ_idx.m1]./new_c.q[circ_idx.m1]
    Dp1[:] = Dp1 - 0.5new_c.α              - W./new_c.q

    D0[:] = D0 + (
                  new_c.r/P.δt
                  + 0.5*(new_c.α - new_c.α[circ_idx.m1])
                  + (W./new_c.q + W[circ_idx.m1]./new_c.q[circ_idx.m1])
                 )

    """
    From the documentation:
    rem(x, y, r::RoundingMode)

      Compute the remainder of x after integer division by y, with the quotient rounded according to the
      rounding mode r. In other words, the quantity

      x - y*round(x/y,r)

      without any intermediate rounding.

        •    if r == RoundNearest, then the result is exact, and in the interval [-|y|/2, |y|/2]. See
            also RoundNearest.

    So that rem.(x, 2π, RoundNearest) gives a float in [-π, π]
    """
    diff_pow_3 = 0.5P.N_kb*(-(rem.(c.θ[circ_idx.p1] - c.θ, 2π, RoundNearest)./c.q).^3 + (rem.(c.θ - c.θ[circ_idx.m1], 2π, RoundNearest)./c.q[circ_idx.m1]).^3)

    f[:] = (new_c.r.*c.θ/P.δt
            + diff_pow_3
            + 0.5*vec(CSC.@dotprod(∇W, (new_c.n[circ_idx.p1,:] + new_c.n)) - CSC.@dotprod(∇W[circ_idx.m1,:], (new_c.n + new_c.n[circ_idx.m1,:])))
           )

    # Periodic boundary conditions:
    f[1] += 2π*(Dm2[1] + Dm1[1])
    f[2] += 2π*Dm2[2]
    f[end]   -= 2π*(Dp2[end] + Dp1[end])
    f[end-1] -= 2π*Dp2[end-1]

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

    finite_differences_3(new_c, Dm2, Dm1, D0, Dp1, Dp2; prefactor=P.N_kb, dual=true)

    Dm1[:] = Dm1 + 1.5P.N_kb*new_c.k.^2 ./new_c.r                           + 0.5new_c.α - W./new_c.r
    Dp1[:] = Dp1 + 1.5P.N_kb*new_c.k[circ_idx.p1].^2 ./new_c.r[circ_idx.p1] - 0.5new_c.α - W./new_c.r[circ_idx.p1]

    D0[:] = new_c.q/P.δt - (Dm2 + Dm1 + Dp1 + Dp2)

    n = new_c.n
    n_p1 = circshift(n, -1)

    area = 0.25*sum(c.q.*(CSC.@dotprod(c.n, c.Y) + CSC.@dotprod(c.n[circ_idx.p1,:], c.Y)))
    PA = -P.N_mu*(area-P.N_target_area)

    f = (
         vec(new_c.q .* c.Y)/P.δt
         -0.5(P.N_P - PA)*vec(new_c.r[circ_idx.p1].*n_p1 + new_c.r.*n)
         -0.5*vec(
                new_c.r[circ_idx.p1] .* vec(CSC.@dotprod(new_c.n[circ_idx.p1,:], ∇W)) .* new_c.n[circ_idx.p1,:]
               +new_c.r .* vec(CSC.@dotprod(new_c.n, ∇W)) .* new_c.n
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

    θ0 = π/P.Nnuc
    θ = collect(π/2 + θ0 .+ (-1:(P.Nnuc-2))*2*θ0)

    nc = NucleusCoords(
        P.N_r_init*[cos.(t) sin.(t)] .+ bc, # Y
        zeros(P.Nnuc), # α
        zeros(P.Nnuc), # β
        1/P.N_r_init * ones(P.Nnuc), # k
        sin(π/P.Nnuc)*2P.N_r_init * ones(P.Nnuc), # r
        sin(π/P.Nnuc)*2P.N_r_init * ones(P.Nnuc), # q
        θ, # θ
        [sin.(θ) -cos.(θ)], # n
        log.(sin(π/P.Nnuc)*2P.N_r_init * ones(P.Nnuc)), # η
        2π*P.N_r_init, # L
        CSC.init_circ_idx(P.Nnuc) # circ_idx
       )

    return nc
end

function update_coords(c::NucleusCoords, new_c::NucleusCoords,
                       potentials::CSC.InteractionPotentials,
                       P::Params, F::Flags, temparrays::CSC.TempArrays6)

    N_W = potentials.N_W
    N_∇W = potentials.N_∇W

    if F.DEBUG
        println("PRE alpha, beta, r, q, K, θ")
        # display([c.α c.β c.r c.q c.k c.θ])
        display([sum(c.α)/P.Nnuc sum(c.β)/P.Nnuc sum(c.r)/P.Nnuc sum(c.q)/P.Nnuc sum(c.k)/P.Nnuc sum(c.θ)/P.Nnuc])
        println()
    end

    update_alphabeta(c, new_c, N_W, N_∇W, P, F)
    update_r(c, new_c, N_W, N_∇W, P, F)
    update_K(c, new_c, N_W, N_∇W, P, F, temparrays)
    update_θ(c, new_c, N_W, N_∇W, P, F, temparrays)
    update_Y(c, new_c, N_W, N_∇W, P, F, temparrays)

    if F.DEBUG
        println("POST alpha, beta, r, q, K, θ")
        # display([new_c.α new_c.β new_c.r new_c.q new_c.k new_c.θ])
        display([sum(new_c.α)/P.Nnuc sum(new_c.β)/P.Nnuc sum(new_c.r)/P.Nnuc sum(new_c.q)/P.Nnuc sum(new_c.k)/P.Nnuc sum(new_c.θ)/P.Nnuc])
        println()
        println("R^j            = ", maximum(new_c.Y[:,1]))
        println("R^(j-1) + δt β = ", maximum(c.Y[:,1])+P.δt*c.β[1])
        println("1/K^j          = ", P.Nnuc/sum(new_c.k))
        println("R(r^(j))       = ", 0.5*(sum(new_c.r)/P.Nnuc)/sin(π/P.Nnuc))
    end
end

end # module Nucleus
