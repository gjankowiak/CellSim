module Wall

using CellSimCommon

export compute_field, compute_walls, check_OOB

"""
The wall is implemented using a penalization method
The corresponding potential is computed as
    f(x, y) = g(x, y) + g(-x, y)
where
    g(x, y) = k(f_width + f_β*pulse(y) - H(x))
    H(x) = -h(x) * log(αx)
and
    pulse is a smooth oscillating function parameterized by
    f_β::Float64      # depth
    f_ω0::Float64     # pulsation
    f_σ::Int          # direction
    f_nk::Int         # number of Fourier components
    (see CellSimCommon.jl)

with this notation,
    H'(x) = -(h'(x) log(αx) + h(x)/x)
    H''(x) = -(h''(x) log(αx) + 2h'(x)/x + h(x)/x^2)

    h is chosen such that h(1/α) = h'(1/α) = h''(1/α) = 0

    a possible such choice is h(x) = min(αx-1, 0)²
"""

macro pulse_init(x)
    return esc(quote
           const σω0 = P.f_σ*P.f_ω0
           pulse = sin.(σω0*$x)
       end)
end

macro dpulse_init(x)
    return esc(quote
           dpulse  =  σω0*cos.(σω0*$x)
           d2pulse = -σω0^2*sin.(σω0*$x)
       end)
end

macro pulse_step(x)
    return esc(quote
           pulse   -= (-1)^i * sin.(σω0*i*$x)/i
       end)
end

macro dpulse_step(x)
    return esc(quote
           dpulse  -=  (-1)^i * σω0  *  cos.(σω0*i*$x)
           d2pulse -= -(-1)^i * σω0^2*i*sin.(σω0*i*$x)
       end)
end

# return esc(:(2/pi*atan.($α*$x))) # atan
# return esc(:(2*$α/pi./(1+($α*$x).^2))) # atan
# return esc(:(-2*($α)^2*$x./(1+($α*$x).^2).^2))
# return esc(:(tan.(pi/2*($y))/$α)) # atan

function g(x::Vector{Float64}, α::Float64)
    return -min.(α*x-1, 0.0).^2 .* log.(α*x)
end

function g_p(x::Vector{Float64}, α::Float64)
    return -(2α*min.(α*x-1, 0.0).*log.(α*x) + min.(α*x-1, 0.0).^2./x)
end

function g_pp(x::Vector{Float64}, α::Float64)
    return -(2α^2*(x.<(1/α)).*log.(α*x) + 4α*min.(α*x-1, 0.0)./x - min.(α*x-1, 0.0).^2./x.^2)
end

macro potential_thres()
    return :(1e-1)
end

function compute_walls(y::Vector, P::Params, F::Flags, levelset::Float64=0.0; right::Bool=false)
    N = size(y)
    @pulse_init(y)
    α = P.f_α
    for i = 2:Int(P.f_nk)
        @pulse_step(y)
    end
    if levelset == 0.0
        if F.circular_wall
            return inv_polar_projection([(1-2*right)*(-P.f_β*pulse-P.f_width)+P.polar_shift y], P)
        else
            return [(1-2*right)*(-P.f_β*pulse-P.f_width)+P.polar_shift y]
        end
    else
        if F.circular_wall
            return inv_polar_projection([(1-2*right)*(-P.f_β*pulse-P.f_width+1/α)+P.polar_shift y], P)
        else
            return [(1-2*right)*(-P.f_β*pulse-P.f_width+1/α)+P.polar_shift y]
        end
    end
end

function compute_field(x::Matrix, P::Params, F::Flags; gradient::Bool=true, hessian::Bool=true, weighted::Bool=false)
    local N = P.N

    virt_x = x

    if F.circular_wall
        (virt_x, DTx, (HT1x, HT2x)) = polar_projection(x, P)
    end

    prefac = 1e0
    α = P.f_α
    β = P.f_β
    if α == 0
        return zeros(N), zeros(N, 2), spzeros(2N, 2N)
    end
    @pulse_init(virt_x[:,2])
    @dpulse_init(virt_x[:,2])
    for i = 2:Int(P.f_nk)
        @pulse_step(virt_x[:,2])
        @dpulse_step(virt_x[:,2])
    end

    x_right  = max.(1e-25, + P.f_width + β*pulse - virt_x[:,1])
    x_left = max.(1e-25, + P.f_width + β*pulse + virt_x[:,1])
    f_right  = g(x_right, α)
    f_left = g(x_left, α)
    local f = f_right + f_left

    if gradient
        ∇f = [g_p(x_left, α)-g_p(x_right, α) β*dpulse.*(g_p(x_right, α)+g_p(x_left, α))]
        if hessian
            DN = β*dpulse.*(g_pp(x_left, α)-g_pp(x_right, α))
            D0 = [g_pp(x_right, α)+g_pp(x_left, α) (β*d2pulse.*(g_p(x_right, α)+g_p(x_left, α))+β^2*dpulse.^2.*(g_pp(x_right, α)+g_pp(x_left, α)))]
            H = spdiagm((DN, D0, DN), (-N, 0, N), 2N, 2N)
            if F.circular_wall
                M∇f1 = spdiagm((∇f[:,1], [∇f[:,1] ∇f[:,1]], ∇f[:,1]), (-N, 0, N), 2N, 2N)
                M∇f2 = spdiagm((∇f[:,2], [∇f[:,2] ∇f[:,2]], ∇f[:,2]), (-N, 0, N), 2N, 2N)
                return f, reshape(DTx'*vec(∇f), N, 2), DTx'*H'*DTx + M∇f1.*HT1x + M∇f2.*HT2x
            else
                return f, ∇f, H
            end
        end
        if F.circular_wall
            return prefac*f, prefac*reshape(DTx'*vec(∇f), N, 2), prefac*spzeros(2N, 2N)
        else
            return prefac*f, prefac*∇f, prefac*spzeros(2N, 2N)
        end
    end
    if weighted
        return prefac*f, Matrix(0, 0), spzeros(2N, 2N)
    else
        return prefac*f, Matrix(0, 0), spzeros(2N, 2N)
    end
end

function check_OOB(x::Matrix, P::Dict{String,Float64})
    N = size(x,1)
    @pulse_init(x[:,2])
    for i = 2_Int(P.f_nk)
        @pulse_step(x[:,2])
    end
    x_right  = - P.f_width - P.f_β*pulse + x[:,1]
    x_left = - P.f_width - P.f_β*pulse - x[:,1]
    return (x_right .< 1e-15) | (x_left .< 1e-15)
end

function polar_projection(x::Matrix, P::Params)
    ymax = 1.0
    N = P.N
    r = sqrt.(sum(abs2, x, 2))
    Tx = [r-P.polar_shift ymax/pi*atan.(x[:,2]./(r+x[:,1]))]
    DTx =  spdiagm((-ymax/(2π)*x[:,2]./r.^2, [x[:,1]./r ymax/(2π)*x[:,1]./r.^2], x[:,2]./r), (-N, 0, N), 2N, 2N)
    HT1x = spdiagm((x[:,1].*x[:,2]./r.^1.5, [x[:,2].^2./r.^1.5 x[:,1].^2./r.^1.5], x[:,1].*x[:,2]./r.^1.5), (-N, 0, N), 2N, 2N)
    HT2x = spdiagm(((x[:,2].^2-x[:,1].^2)./r.^4, [-x[:,1].*x[:,2]./r.^4 x[:,1].*x[:,2]./r.^4], (x[:,2].^2-x[:,1].^2)./r.^4), (-N, 0, N), 2N, 2N).*ymax/π
    return (Tx, DTx, (HT1x, HT2x))
end

function inv_polar_projection(x::Matrix, P::Params)
    ymax = 1.0
    return [x[:,1].*cospi.(2x[:,2]/ymax) x[:,1].*sinpi.(2x[:,2]/ymax)]
end


end # module
