module CellSimCommon

export Params, Flags

import SparseArrays
const SA = SparseArrays

function dump_struct(s, filename)
    dump = open(filename, "w")
    for f in fieldnames(typeof(s))
        write(dump, string(f, "\n"))
        show(dump, "text/plain", getfield(s, f))
        write(dump, "\n\n")
    end
    close(dump)
end

struct CircIdx
    m2::Array{Int64,1}
    m1::Array{Int64,1}
    p1::Array{Int64,1}
    p2::Array{Int64,1}
end

function init_circ_idx()
    CircIdx([], [], [], [])
end

function init_circ_idx(N::Int64)
    CircIdx(circshift(1:N,  2), circshift(1:N,  1),
             circshift(1:N, -1), circshift(1:N, -2))
end

struct Params
    M::Int64      # max. number of iterations

    Δσ::Float64   # space step
    δt::Float64   # time step

    # Cortex related parameters
    N::Int64      # number of points on the cortex
    P::Float64    # pressure
    K::Float64    # membrane elasticity
    Ka::Float64   # cortex viscosity
    c::Float64    # polymerization speed

    # Initial condition parameters
    x0_a::Float64      # initial ellipsis width
    x0_b::Float64      # initial ellipsis height
    x0_shift::Float64  # initial vertical shift

    # Confinement field
    f_α::Float64      # sharpness
    f_β::Float64      # depth
    f_ω0::Float64     # pulsation
    f_σ::Int          # direction
    f_nk::Int         # number of Fourier components
                      # to approximate a saw-tooth signal
    f_width::Float64  # mean width
    f_iwidth::Float64 # inner width

    drag_gauss_power::Float64 # gaussian power in the drag mask
    drag_gauss_width::Float64 # gaussian width in the drag mask

    mass_gauss_power::Float64 # gaussian power in the mass mask
    mass_gauss_width::Float64 # gaussian width in the mass mask

    polar_shift::Float64 # how much to shift the obstacle to the right before polar projection

    # Centrosome related parameters
    k_MT::Float64 # (isotropic) MT/cortex friction coefficient

    MT_potential_power::Float64 # α, where the MT force is F_MT = k (Xi-Xc)*|Xi-Xc|^α
    MT_factor::Float64 # the prefactor k in F_MT

    # Nucleus related parameters
    Nnuc::Int64   # number of points on the nucleus
    N_P::Float64  # pressure
    N_kb::Float64 # bending stiffness
    N_ω::Float64  # inplane stiffness
    N_kc::Float64 # centrosome link stiffness
    N_l0c::Float64 # centrosome link rest length
    N_r_init::Float64 # initial nucleus radius
end

struct Flags
    # model options
    confine::Bool
    adjust_drag::Bool
    polymerize::Bool
    continuous::Bool
    circular_wall::Bool
    cortex::Bool
    centrosome::Bool
    nucleus::Bool
    weighted_confinement::Bool

    # scheme/solver options
    innerloop::Bool

    # visualization options
    plot::Bool
    pretty::Bool
    landscape_plot::Bool
    follow_cam::Bool
    plot_drag::Bool

    # output options
    dryrun::Bool
    write_animation::Bool
end

struct InteractionPotentials
    N_W::Vector{Float64}
    N_∇W::Matrix{Float64}
    C_∇W::Matrix{Float64}
    CS_∇W::Matrix{Float64}
end

struct TempArrays6
    v1::Vector{Float64}
    v2::Vector{Float64}
    v3::Vector{Float64}
    v4::Vector{Float64}
    v5::Vector{Float64}
    v6::Vector{Float64}
end

struct TempArrays12
    v1::Vector{Float64}
    v2::Vector{Float64}
    v3::Vector{Float64}
    v4::Vector{Float64}
    v5::Vector{Float64}
    v6::Vector{Float64}
    v7::Vector{Float64}
    v8::Vector{Float64}
    v9::Vector{Float64}
    v10::Vector{Float64}
    v11::Vector{Float64}
    v12::Vector{Float64}
end

macro ta6_tuple(ta)
    return esc(:(($ta.v1, $ta.v2, $ta.v3, $ta.v4, $ta.v5, $ta.v6)))
end

macro ta12_tuple(ta)
    return esc(:(($ta.v1, $ta.v2, $ta.v3, $ta.v4, $ta.v5, $ta.v6,
                  $ta.v7, $ta.v8, $ta.v9, $ta.v10, $ta.v11, $ta.v12)))
end

struct Plotables
    field::Vector{Float64}
    ∇field::Matrix{Float64}
    transport_force::Matrix{Float64}
    drag_force::Matrix{Float64}
    mass_source::Vector{Float64}
    mass_source_int::Vector{Float64}
    mt_force::Vector{Float64}
    mt_force_indiv::Matrix{Float64}
end

function new_plotables(N::Int64)
    return Plotables(
        zeros(N),  # field
        zeros(N,2), # ∇field
        zeros(N,2), # transport_force
        zeros(N,2), # drag_force
        zeros(N),   # mass_source
        zeros(N),   # mass_source_int
        zeros(2),   # mt_force
        zeros(N+1,2),# mt_force_indiv
       )
end

struct Metrics
    iter::Int64
end

macro looped(array, i_min, i_max)
    return :( $array[mod((($i_min)-1):(($i_max)-1), size($array, 1))+1] )
end

macro delta(x)
    return :(circshift($x, -1) - $x)
end

macro delta!(src, dst)
    return esc(:($dst[:] = circshift!($dst, $src, -1) - $src))
end

macro entry_norm(x)
    return esc(:(vec(sqrt.(sum(abs2, $x; dims=2)))))
end

macro barycenter(x)
    return :(sum(view($x,:,2))/size($x,1))
end

macro dotprod(x,y)
    return esc(:(sum($x.*$y; dims=2)))
end

macro new_point(x_max, x_nb)
    return :(($x_max+$x_nb)/2)
end

macro det(x1, x2)
    return :($x1[:,1].*$x2[:,2] - $x1[:,2].*$x2[:,1])
end

macro repdiag(M, n)
    return :( SA.spdiagm(0 => repeat($M, $n, 1)) )
end

macro repdiagblk(M, n)
    return :( SA.blockdiag(ntuple((_) -> $M, $n)...) )
end

@inbounds function shift!(src::Array{Float64}, dst::Array{Float64})
    last = src[end,:]
    dst[2:end,:] = src[1:end-1,:]
    dst[1,:] = last
    return dst
end

@inbounds function ishift!(src::Array{Float64}, dst::Array{Float64}, offset::Int=1)
    first = src[1:offset,:]
    dst[1:end-offset,:] = src[1+offset:end,:]
    dst[end+1-offset:end,:] = first
    return dst
end

function perp!(src::Matrix, dst::Matrix)
    if src === dst
        throw("Source and destination arrays must be different")
    end
    dst[:,1] = src[:,2]
    dst[:,2] = -src[:,1]
end

macro bc_scalar(v)
    return esc(:(repeat($v, 2, 1)))
end

"""
    pointwise_projection(v)

    returns the matrix corresponding to the projection on v: x → (x∙v) v
"""
function pointwise_projection(v)
    N = size(v, 1)
    return SA.spdiagm( 0 => [v[:,1].*v[:,1]; v[:,2].*v[:,2]],
                       -N => v[:,1].*v[:,2],
                        N => v[:,1].*v[:,2])
    # return SA.spdiagm(([v[:,1].*v[:,1]; v[:,2].*v[:,2]], v[:,1].*v[:,2], v[:,1].*v[:,2]), (0, -N, N), 2*N, 2*N)
end

"""
    pointwise_projection(bra, ket)

    returns the matrix corresponding to the operator x → (x∙bra) ket
"""
function pointwise_projection(bra, ket)
    N = size(bra, 1)
    return SA.spdiagm( 0 => [bra[:,1].*ket[:,1]; bra[:,2].*ket[:,2]],
                       -N => bra[:,1].*ket[:,2],
                        N => ket[:,1].*bra[:,2])
    # return SA.spdiagm(([bra[:,1].*ket[:,1]; bra[:,2].*ket[:,2]], bra[:,1].*ket[:,2], ket[:,1].*bra[:,2]), (0, -N, N), 2*N, 2*N)
end

"""
    pointwise_dot_prod(v)

    returns the matrix corresponding to the operator x → (x∙v)
"""
function pointwise_dot_prod(v)
    N = size(v, 1)
    return [v[:, 1].*SA.sparse(SA.I, N, N) v[:,2].*SA.sparse(SA.I, N, N)]
end

"""
Stolen from https://github.com/mlubin/NaNMath.jl/blob/master/src/NaNMath.jl

NaNMath.sum(A)
##### Args:
* `A`: An array of floating point numbers
##### Returns:
*    Returns the sum of all elements in the array, ignoring NaN's.
##### Examples:
```julia
using NaNMath as nm
nm.sum([1., 2., NaN]) # result: 3.0
```
"""
function nansum(x::AbstractArray{T}) where T<:AbstractFloat
    if length(x) == 0
        result = zero(eltype(x))
    else
        result = convert(eltype(x), NaN)
        for i in x
            if !isnan(i)
                if isnan(result)
                    result = i
                else
                    result += i
                end
            end
        end
    end

    if isnan(result)
        Base.warn_once("All elements of the array, passed to \"sum\" are NaN!")
    end
    return result
end

end # module
