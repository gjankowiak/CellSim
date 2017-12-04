module Wall

using AdhCommon

export compute_field, compute_walls, check_OOB

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

macro potential_thres()
    return :(5e-3)
end

macro barrier_p()
    return :(1)
end

function barrier_f(x)
    return log.(x)
end

function barrier_fp(x)
    return 1./x
end

function barrier_fpp(x)
    return -1./x.^2
end

function barrier_full(x)
    const p = @barrier_p
    if p == 1
        return barrier_f(x)
    else
        return barrier_f(x) .* abs.(barrier_f(x)).^(p-1)
    end
end

function inv_barrier_full(x)
    const p = @barrier_p
    if p == 1
        return exp(-x)
    else
        return exp(-x.^(1/p))
    end
end

function barrier_full_p(x)
    const p = @barrier_p
    if p == 1
        return barrier_fp(x)
    else
        f = barrier_f(x)
        return p * barrier_fp(x) .* abs.(f).^(p-1)
    end
end

function barrier_full_pp(x)
    const p = @barrier_p
    if p == 1
        return barrier_fpp(x)
    else
        f = barrier_f(x)
        return p * abs.(f).^(p - 2) .* (barrier_fpp(x) .* abs.(f) + (p-1)*sign(f).*barrier_fp(x).^2)
    end
end

function compute_walls(y::Vector, params::Params, levelset::Float64=0.0)
    N = size(y)
    @pulse_init(y)
    threshold_f = @potential_thres
    width_active_wall = 3e-2*params.f_width
    α = -threshold_f/barrier_full([width_active_wall])[1]
    α = params.f_α
    for i = 2:Int(params.f_nk)
        @pulse_step(y)
    end
    if levelset == 0.0
        return params.f_β*pulse + params.f_width
    else
        return inv_barrier_full(threshold_f./α)-params.f_β*pulse - params.f_width
    end
end

function compute_field(x::Matrix, params::Params; gradient::Bool=true, hessian::Bool=true)
    local N = params.N
    threshold_f = @potential_thres
    width_active_wall = 3e-2*params.f_width
    α = -threshold_f/barrier_full(width_active_wall)
    α = params.f_α
    if α == 0
        return zeros(N), zeros(N, 2), spzeros(2N, 2N)
    end
    @pulse_init(x[:,2])
    @dpulse_init(x[:,2])
    for i = 2:Int(params.f_nk)
        @pulse_step(x[:,2])
        @dpulse_step(x[:,2])
    end
    x_left  = max.(1e-25, params.f_width + params.f_β*pulse + x[:,1])
    x_right = max.(1e-25, params.f_width + params.f_β*pulse - x[:,1])
    f_left  = -α*barrier_full(x_left)
    f_right = -α*barrier_full(x_right)
    local f = f_left + f_right

    trim_mask = (f_left .> -α*barrier_full(width_active_wall)) .| (f_right .> -α*barrier_full(width_active_wall))

    if gradient
        ∇f = -α*trim_mask.*[barrier_full_p(x_left)-barrier_full_p(x_right) params.f_β*dpulse.*(barrier_full_p(x_left)+barrier_full_p(x_right))]
        if hessian
            DN = params.f_β*dpulse.*trim_mask.*(barrier_full_pp(x_left)-barrier_full_pp(x_right))
            D0 = [trim_mask.*(barrier_full_pp(x_left)+barrier_full_pp(x_right)) trim_mask.*(params.f_β*d2pulse.*(barrier_full_p(x_left)+barrier_full_p(x_right))+params.f_β^2*dpulse.^2.*(barrier_full_pp(x_left)+barrier_full_pp(x_right)))]
            H = -α*spdiagm((DN, D0, DN), (-N, 0, N), 2N, 2N)
            return f, ∇f, H
        end
        return f, ∇f, spzeros(2N, 2N)
    end
    return f, Matrix(0, 0), spzeros(2N, 2N)
end

function check_OOB(x::Matrix, params::Dict{String,Float64})
    N = size(x,1)
    @pulse_init(x[:,2])
    for i = 2_Int(params.f_nk)
        @pulse_step(x[:,2])
    end
    x_left  = params.f_width + params.f_β*pulse + x[:,1]
    x_right = params.f_width + params.f_β*pulse - x[:,1]
    return (x_left .< 1e-15) | (x_right .< 1e-15)
end


end # module
