
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

# max_strain_criteria: consumes local engineering strains [ε1, ε2, γ12] returned by
# local_deformations (third component is engineering shear strain γ12 = 2ε12).
# Allowable strains are derived from stress allowables and elastic constants:
#   ε1t = Xt/E1,  ε1c = Xc/E1
#   ε2t = Yt/E2,  ε2c = Yc/E2
#   γ12_allow = S12/G12
function max_strain_criteria(lam::Laminate, load::AbstractVector{<:Real})
    strain = local_deformations(lam, load)

    crit = Vector{Vector{Bool}}()

    for i in eachindex(lam.sheets)
        m = lam.sheets[i].material
        check_material_failure_prop(m)
        ϵ = strain[i]
        ε1t = m.Xt / m.E1
        ε1c = m.Xc / m.E1
        ε2t = m.Yt / m.E2
        ε2c = m.Yc / m.E2
        γ12_allow = m.S12 / m.G12
        push!(crit, [-ε1c < ϵ[1] < ε1t, -ε2c < ϵ[2] < ε2t, abs(ϵ[3]) < γ12_allow])
    end

    return crit
end

function tsai_hill_criteria(lam::Laminate, load::AbstractVector{<:Real})
    tensions = local_tensions(lam, load)

    crit = Float64[]

    for i in eachindex(lam.sheets)
        m = lam.sheets[i].material
        check_material_failure_prop(m)
        σ = tensions[i]
        σ1, σ2, σ3 = σ
        if σ1 > 0 && σ2 > 0
            FI = (σ1 / m.Xt)^2 - (σ1 / m.Xt) * (σ2 / m.Xt) + (σ2 / m.Yt)^2 + (σ3 / m.S12)^2
        elseif σ1 < 0 && σ2 > 0
            FI = (σ1 / m.Xc)^2 + (σ1 / m.Xc) * (σ2 / m.Xc) + (σ2 / m.Yt)^2 + (σ3 / m.S12)^2
        elseif σ1 > 0 && σ2 < 0
            FI = (σ1 / m.Xt)^2 + (σ1 / m.Xt) * (σ2 / m.Xt) + (σ2 / m.Yc)^2 + (σ3 / m.S12)^2
        else # σ1 <= 0 && σ2 <= 0 (or both zero)
            FI = (σ1 / m.Xc)^2 - (σ1 / m.Xc) * (σ2 / m.Xc) + (σ2 / m.Yc)^2 + (σ3 / m.S12)^2
        end
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
