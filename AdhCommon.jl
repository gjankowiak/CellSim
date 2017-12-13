module AdhCommon

export Params, Flags

function init(N::Int64)
    const global circ_idx_m2 = circshift(1:N,  2)
    const global circ_idx_m1 = circshift(1:N,  1)
    const global circ_idx_p1 = circshift(1:N, -1)
    const global circ_idx_p2 = circshift(1:N, -2)
end

immutable Params
    N::Int64      # number of points
    M::Int64      # max. number of iterations

    Δσ::Float64   # space step
    δt::Float64   # time step

    P::Float64    # pressure
    K::Float64    # membrane elasticity
    Ka::Float64   # cortex viscosity
    c::Float64    # polymerization speed

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
end

immutable Flags
    confine::Bool
    adjust_drag::Bool
    polymerize::Bool
    dryrun::Bool
    plot::Bool
    pretty::Bool
    continuous::Bool
    innerloop::Bool
    weighted_confinement::Bool
end

immutable Metrics
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
    return esc(:(vec(sqrt.(sum(abs2, $x, 2)))))
end

macro barycenter(x)
    return :(sum(view($x,:,2))/size($x,1))
end

macro dotprod(x,y)
    return esc(:(sum($x.*$y,2)))
end

macro new_point(x_max, x_nb)
    return :(($x_max+$x_nb)/2)
end

macro det(x1, x2)
    return :($x1[:,1].*$x2[:,2] - $x1[:,2].*$x2[:,1])
end

macro repdiag(M, n)
    return :( spdiagm(repmat($M, $n, 1), 0) )
end

macro repdiagblk(M, n)
    return :( blkdiag(ntuple((_) -> $M, $n)...) )
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
    if src === b
        throw("Source and destination arrays must be different")
    end
    dst[:,1] = src[:,2]
    dst[:,2] = -src[:,1]
end

macro bc_scalar(v)
    return esc(:(repmat($v, 2, 1)))
end

"""
    pointwise_projection(v)

    returns the matrix corresponding to the projection on v: x → (x∙v) v
"""
function pointwise_projection(v)
    N = size(v, 1)
    return spdiagm(([v[:,1].*v[:,1]; v[:,2].*v[:,2]], v[:,1].*v[:,2], v[:,1].*v[:,2]), (0, -N, N), 2*N, 2*N)
end

"""
    pointwise_projection(bra, ket)

    returns the matrix corresponding to the operator x → (x∙bra) ket
"""
function pointwise_projection(bra, ket)
    N = size(bra, 1)
    return spdiagm(([bra[:,1].*ket[:,1]; bra[:,2].*ket[:,2]], bra[:,1].*ket[:,2], ket[:,1].*bra[:,2]), (0, -N, N), 2*N, 2*N)
end

end
