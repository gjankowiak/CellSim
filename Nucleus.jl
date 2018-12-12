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

    # local length element
    η::Vector{Float64}

    # nuclear membrane lenght
    L::Float64

    circ_idx::CircIdx
end

function compute_contact_force(cor_coords::PointCoords, nuc_coords::NucleusCoords,
                               P::Params, F::Flags)
end

function compute_centronuclear_force(cor_coords::PointCoords, nuc_coords::NucleusCoords,
                                     P::Params, F::Flags)
end

function update_alphabeta(nuc_coords::NucleusCoords, new_nuc_coords::NucleusCoords,
                          W::Vector{Float64}, ∇W::Matrix{Float64},
                          P::Params, F::Flags)
    B = sum(nuc_coords.r .* nuc_coords.k .* nuc_coords.β) / nuc_coords.L
    new_nuc_coords.α[1] = 0
    for i in 2:P.Nnuc
        new_nuc_coords.α[i] = (new_nuc_coords.α[i-1]
                           + nuc_coords.r[i]*(-nuc_coords.k[i]*nuc_coords.β[i]
                                              +nuc_coords.B)
                           + (nuc_coords.L/P.Nnuc - nuc_coords.r[i])*P.N_ω)
    end
end

function update_r(nuc_coords::NucleusCoords, new_nuc_coords::NucleursCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags)
    new_nuc_coords.η[:] = nuc_coords.η + P.δt*(nuc_coords.k .* nuc_coords.β
                                               +(new_nuc_coords.α - new_nuc_coords.α[new_nuc_coords.circ_idx.m1])
                                               ./nuc_coords.r)
    new_nuc_coords.r[:] = exp.(new_nuc_coords.η)
end

function update_K(nuc_coords::NucleusCoords, new_nuc_coords::NucleursCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags)
end

function update_υ(nuc_coords::NucleusCoords, new_nuc_coords::NucleursCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags)
end

function update_Y(nuc_coords::NucleusCoords, new_nuc_coords::NucleursCoords,
                  W::Vector{Float64}, ∇W::Matrix{Float64},
                  P::Params, F::Flags)
end


end # module Nucleus
