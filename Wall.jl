module Wall

using AdhCommon

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
    (see AdhCommon.jl)

with this notation,
    H'(x) = -(h'(x) log(αx) + h(x)/x)
    H''(x) = -(h''(x) log(αx) + 2h'(x)/x + h(x)/x^2)

    h is chosen such that h(1/α) = h'(1/α) = h''(1/α) = 0

    a possible such choice is h(x) = min(αx-1, 0)²
"""

macro pulse_init(x)
    return esc(quote
           const σω0 = params.f_σ*params.f_ω0
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
    return -(2α^2*(x.<(1/α)).*log.(α*x) + 4α*min.(α*x-1, 0.0)./x - min.(α*x-1, 0.0).^2/x.^2)
end

macro potential_thres()
    return :(1e-1)
end

function compute_walls(y::Vector, params::Params, levelset::Float64=0.0)
    N = size(y)
    @pulse_init(y)
    α = params.f_α
    for i = 2:Int(params.f_nk)
        @pulse_step(y)
    end
    if levelset == 0.0
        return -params.f_β*pulse - params.f_width
    else
        return -params.f_β*pulse - params.f_width + 1/α
    end
end

function compute_field(x::Matrix, params::Params; gradient::Bool=true, hessian::Bool=true)
    local N = params.N
    prefac = 1e1
    α = params.f_α
    β = params.f_β
    if α == 0
        return zeros(N), zeros(N, 2), spzeros(2N, 2N)
    end
    @pulse_init(x[:,2])
    @dpulse_init(x[:,2])
    for i = 2:Int(params.f_nk)
        @pulse_step(x[:,2])
        @dpulse_step(x[:,2])
    end
    x_right  = max.(1e-25, + params.f_width + β*pulse - x[:,1])
    x_left = max.(1e-25, + params.f_width + β*pulse + x[:,1])
    f_right  = g(x_right, α)
    f_left = g(x_left, α)
    local f = f_right + f_left

    if gradient
        ∇f = [g_p(x_left, α)-g_p(x_right, α) β*dpulse.*(g_p(x_right, α)+g_p(x_left, α))]
        if hessian
            DN = β*dpulse.*(g_pp(x_left)-g_pp(x_right))
            D0 = [g_pp(x_right, α)+g_pp(x_left, α) (β*d2pulse.*(g_p(x_right, α)+g_p(x_left, α))+β^2*dpulse.^2.*(g_pp(x_right, α)+g_pp(x_left, α)))]
            H = spdiagm((DN, D0, DN), (-N, 0, N), 2N, 2N)
            return f, ∇f, H
        end
        return prefac*f, prefac*∇f, prefac*spzeros(2N, 2N)
    end
    return prefac*f, Matrix(0, 0), spzeros(2N, 2N)
end

function check_OOB(x::Matrix, params::Dict{String,Float64})
    N = size(x,1)
    @pulse_init(x[:,2])
    for i = 2_Int(params.f_nk)
        @pulse_step(x[:,2])
    end
    x_right  = - params.f_width - params.f_β*pulse + x[:,1]
    x_left = - params.f_width - params.f_β*pulse - x[:,1]
    return (x_right .< 1e-15) | (x_left .< 1e-15)
end


end # module
