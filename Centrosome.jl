module Centrosome

function compute_angles(X::Array{Float64,2})
    return angle.(X[:,1] + X[:,2]*im)
end

function intercept(α::Array{Float64,1}, α_s::SubArray{Float64,1}, α_min::Array{Float64,1}, α_max::Array{Float64,1}, α_Ox::BitArray{1}, θ::Float64)
    return find(((α .<= θ .<= α_s) .& !α_Ox) .| (α_Ox .& ((θ .< α_min) .| (θ .> α_max))))
end

function closest_intercept(X::Array{Float64,2}, X_s::SubArray{Float64,2}, icpt::Array{Int64,1}, θ::Float64)

    # helper matrices
    M = [-sin(θ); cos(θ)]
    MM = [0 1; -1 0]

    # compute the distance between the center (here 0) and the intersection point
    μ = sum(X[icpt,:] * MM .* X_s[icpt,:], 2)./((X_s[icpt,:] - X[icpt,:])*M)

    # take the smallest among all intersections
    r, _min_idx = findmin(μ)

    # get the corresponding index in X
    min_idx = icpt[_min_idx]

    # compute the wieght of the convex combination
    # avoid small denominators first
    if abs(X[min_idx,1] - X_s[min_idx,1]) < 1e-8
        dim = 2
        trig = sin(θ)
    else
        dim = 1
        trig = cos(θ)
    end

    λ = (r*trig-X[min_idx,dim])./(X_s[min_idx,dim]-X[min_idx,dim])
    return min_idx, λ, r
end

function compute_intercept(X::Array{Float64,2}, center::Array{Float64,2}, angles::Array{Float64,1})
    N = size(X, 1)

    # do all computations with the center of radiation at the origin
    Xc = X.-center

    # compute the apparent angle of each point relative to the center
    α = compute_angles(Xc)

    # shifted vectors
    idx_s = circshift(1:N, -1)
    Xc_s = view(Xc, idx_s, :)
    α_s = view(α, idx_s)
    α_min, α_max = min.(α, α_s), max.(α, α_s)
    α_Ox = abs(α-α_s) .> π

    # storage
    # the first point of the interected segment at angle angles[i]
    idx_Xθ = zeros(angles, Int64)

    # the weight of the convex combination between X[idx_Xθ[i]] and
    # X[idx_Xθ[i]+1] which defines the intersection point with the
    # ray of angle angles[i]
    λ = zeros(angles)

    # the distance of the intersection point with the center
    μ = zeros(angles)

    for (i, θ) in enumerate(angles)
        # compute all intersected segments at angle θ
        icpt = intercept(α, α_s, α_min, α_max, α_Ox, θ)
        # find the intersection point closest to the center
        # and corresponding weights and radius
        idx_Xθ[i], λ[i], μ[i] = closest_intercept(Xc, Xc_s, icpt, θ)
    end

    return idx_Xθ, λ, μ
end

end # module
