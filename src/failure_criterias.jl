
function max_tensions_criteria(lam::Laminate, load::AbstractVector{<:Real})
    tensions = local_tensions(lam, load)

    crit = Vector{Vector{Bool}}()

    for i in eachindex(lam.sheets)
        m = lam.sheets[i].material
        check_material_failure_prop(m)
        σ = tensions[i]
        push!(crit, [-m.Xc < σ[1] < m.Xt, -m.Yc < σ[2] < m.Yt, abs(σ[3]) < m.S12])
    end

    return crit
end

function max_strain_criteria(lam::Laminate, load::AbstractVector{<:Real})
    strain = local_deformations(lam, load)

    crit = Vector{Vector{Bool}}()

    for i in eachindex(lam.sheets)
        m = lam.sheets[i].material
        check_material_failure_prop(m)
        ϵ = strain[i]
        push!(crit, [-m.Xc < ϵ[1] < m.Xt, -m.Yc < ϵ[2] < m.Yt, abs(ϵ[3]) < m.S12])
    end

    return crit
end

"""
    tsai_hill_criteria(lam::Laminate, load::AbstractVector{<:Real})

Compute the Tsai–Hill failure index for each ply using local ply stresses
`(σ1, σ2, τ12)` and sign-dependent allowables:

`FI = (σ1/X)^2 - σ1σ2/(X*Y) + (σ2/Y)^2 + (τ12/S12)^2`

where `X = Xt` for `σ1 ≥ 0` and `X = Xc` for `σ1 < 0`, while
`Y = Yt` for `σ2 ≥ 0` and `Y = Yc` for `σ2 < 0`.
This applies the same equation in all four `(σ1, σ2)` sign quadrants,
with tensile/compressive allowables selected by stress sign.
"""
function tsai_hill_criteria(lam::Laminate, load::AbstractVector{<:Real})
    tensions = local_tensions(lam, load)

    crit = Float64[]

    for i in eachindex(lam.sheets)
        m = lam.sheets[i].material
        check_material_failure_prop(m)
        σ = tensions[i]
        σ1, σ2, σ3 = σ
        X = σ1 < 0 ? m.Xc : m.Xt
        Y = σ2 < 0 ? m.Yc : m.Yt
        FI = (σ1 / X)^2 - (σ1 * σ2) / (X * Y) + (σ2 / Y)^2 + (σ3 / m.S12)^2
        push!(crit, FI)
    end

    return crit
end

function tsai_wu_criteria(lam::Laminate, load::AbstractVector{<:Real})
    tensions = local_tensions(lam, load)

    crit = Float64[]

    for i in eachindex(lam.sheets)
        m = lam.sheets[i].material

        check_material_failure_prop(m)

        F1 = 1 / m.Xt - 1 / m.Xc
        F2 = 1 / m.Yt - 1 / m.Yc
        F11 = 1 / (m.Xt * m.Xc)
        F22 = 1 / (m.Yt * m.Yc)
        F66 = 1 / (m.S12^2)
        F12 = -sqrt(F11 * F22) / 2

        σ = tensions[i]
        σ1, σ2, σ6 = σ

        FI =
            (F1 * σ1 + F2 * σ2) + (F11 * σ1^2 + F22 * σ2^2 + F66 * σ6^2) + 2 * F12 * σ1 * σ2

        push!(crit, FI)
    end

    return crit
end

function check_material_failure_prop(m::OrthotropicMaterial)
    if isnothing(m.Xt) ||
        isnothing(m.Yt) ||
        isnothing(m.Xc) ||
        isnothing(m.Yc) ||
        isnothing(m.S12)
        error("Some Material properties are undefined.")
    end

    return true
end
