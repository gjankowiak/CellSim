module Masks

using AdhCommon

function compute_drag_mask(x::Array, half_w::Int, pow::Float64=2.0)
    const N = size(x, 1)

    const x_min, x_min_idx = findmin(view(x,:,2))
    const x_max, x_max_idx = findmax(view(x,:,2))

    mask = zeros(N)
    profile = exp(-linspace(-3, 3, 2*half_w+1).^pow)
    profile .*= 0.5/sum(profile)

    AdhCommon.@looped(mask, x_min_idx-half_w, x_min_idx+half_w) = profile
    AdhCommon.@looped(mask, x_max_idx-half_w, x_max_idx+half_w) = profile

    return mask
end

function compute_drag_mask(x::Array, ΔL::Array, half_w::Float64, pow::Float64=2.0)
    const N = size(x, 1)

    const x_min, x_min_idx = findmin(view(x,:,2))
    const x_max, x_max_idx = findmax(view(x,:,2))

    dst_from_min = split_cumsum(ΔL, x_min_idx)
    dst_from_max = split_cumsum(ΔL, x_max_idx)

    profile_min = exp(-(dst_from_min/half_w).^pow)
    profile_min .*= profile_min .> 1e-5
    profile_min .*= 0.5/sum(profile_min)
    profile_max = exp(-(dst_from_max/half_w).^pow)
    profile_max .*= profile_max .> 1e-5
    profile_max .*= 0.5/sum(profile_max)

    return profile_min + profile_max
end

#function compute_mask(N::Int)
    #@polymerization_bounds

    #mask = ones(N)
    #const mask_width = round(Int, 50*N/400)
    #mask[i_start-mask_width:i_start+mask_width] = [abs(i)/mask_width for i in -mask_width:mask_width]
    #mask[i_end-mask_width:i_end+mask_width] = mask[i_start-mask_width:i_start+mask_width]
    #return mask
#end

end # module
