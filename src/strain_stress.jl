
function global_deformations(lam::Laminate, load::AbstractVector{<:Real})
    deform = ABD(lam) \ load
    return deform[1:3], deform[4:6]
end

function global_deformations_per_sheet(lam::Laminate, load::AbstractVector{<:Real})
    ϵ, κ = global_deformations(lam, load)
    t_k = t_pos(lam)
    deforms = Vector{Vector{Float64}}()
    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        push!(deforms, ϵ .+ (t_k[i] + sh.t_i / 2) .* κ)
    end
    return deforms
end

function global_tensions_per_sheet(lam::Laminate, load::AbstractVector{<:Real})
    g_deforms = global_deformations_per_sheet(lam, load)

    tensions = Vector{Vector{Float64}}()

    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        push!(tensions, Qxy(sh) * g_deforms[i])
    end

    return tensions
end

function local_deformations(lam::Laminate, load::AbstractVector{<:Real})
    g_deforms = global_deformations_per_sheet(lam, load)

    deforms = Vector{Vector{Float64}}()

    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        push!(deforms, Tϵ(sh.θ_i) * g_deforms[i])
    end

    return deforms
end

function local_tensions(lam::Laminate, load::AbstractVector{<:Real})
    g_tensions = global_tensions_per_sheet(lam, load)

    tensions = Vector{Vector{Float64}}()

    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        push!(tensions, Tσ(sh.θ_i) * g_tensions[i])
    end

    return tensions
end
