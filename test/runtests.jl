
using Laminates
using Laminates:
    OrthotropicMaterial,
    Laminate,
    Sheet,
    cross_ply,
    angle_ply,
    thickness,
    Q12,
    S12,
    Qxy,
    A,
    B,
    D,
    ABD,
    abd,
    global_deformations,
    global_deformations_per_sheet,
    global_tensions_per_sheet,
    local_deformations,
    local_tensions,
    max_tensions_criteria,
    max_strain_criteria,
    tsai_hill_criteria,
    tsai_wu_criteria,
    Tσ,
    Tϵ
using Test
using Aqua
using JuliaFormatter

function z_interfaces(ts::AbstractVector{<:Real})
    h = sum(ts)
    z = Float64[-h / 2]
    for t in ts
        push!(z, z[end] + t)
    end
    return z
end

identity_matrix(n::Integer) = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]

function q12_reference(mat::OrthotropicMaterial)
    ν21 = mat.ν12 * mat.E2 / mat.E1
    den = 1 - mat.ν12 * ν21
    q11 = mat.E1 / den
    q22 = mat.E2 / den
    q12 = mat.ν12 * mat.E2 / den
    q66 = mat.G12
    return [
        q11 q12 0.0
        q12 q22 0.0
        0.0 0.0 q66
    ]
end

function s12_reference(mat::OrthotropicMaterial)
    return [
        1 / mat.E1 -mat.ν12 / mat.E1 0.0
        -mat.ν12 / mat.E1 1 / mat.E2 0.0
        0.0 0.0 1 / mat.G12
    ]
end

function qbar_reference(mat::OrthotropicMaterial, θ::Real)
    qref = q12_reference(mat)
    q11, q12, q22, q66 = qref[1, 1], qref[1, 2], qref[2, 2], qref[3, 3]
    θr = deg2rad(θ)
    m, n = cos(θr), sin(θr)
    m2, n2 = m^2, n^2
    m3, n3 = m^3, n^3
    m4, n4 = m^4, n^4
    q̄11 = q11 * m4 + 2 * (q12 + 2 * q66) * m2 * n2 + q22 * n4
    q̄22 = q11 * n4 + 2 * (q12 + 2 * q66) * m2 * n2 + q22 * m4
    q̄12 = (q11 + q22 - 4 * q66) * m2 * n2 + q12 * (m4 + n4)
    q̄16 = (q11 - q12 - 2 * q66) * m3 * n - (q22 - q12 - 2 * q66) * m * n3
    q̄26 = (q11 - q12 - 2 * q66) * m * n3 - (q22 - q12 - 2 * q66) * m3 * n
    q̄66 = (q11 + q22 - 2 * q12 - 2 * q66) * m2 * n2 + q66 * (m4 + n4)
    return [
        q̄11 q̄12 q̄16
        q̄12 q̄22 q̄26
        q̄16 q̄26 q̄66
    ]
end

function abd_reference(
    θs::AbstractVector{<:Real}, ts::AbstractVector{<:Real}, mat::OrthotropicMaterial
)
    z = z_interfaces(ts)
    a = zeros(3, 3)
    b = zeros(3, 3)
    d = zeros(3, 3)
    for i in eachindex(θs)
        qbar = qbar_reference(mat, θs[i])
        z0, z1 = z[i], z[i + 1]
        a .+= qbar * (z1 - z0)
        b .+= qbar * (z1^2 - z0^2) / 2
        d .+= qbar * (z1^3 - z0^3) / 3
    end
    return a, b, d
end

@testset "Laminates.jl" begin
    mat_cfrp = OrthotropicMaterial(70e9, 20e9, 0.3, 10e8)
    mat_cfrp_limit = OrthotropicMaterial(70e9, 20e9, 0.3, 10e8, 100, 100, 100, 100, 20)

    lam_cross = cross_ply(10, 0.1, mat_cfrp)

    @test length(lam_cross) == 10
    @test isapprox(Laminates.thickness(lam_cross), 1.0, atol=1e-15)

    lam_angle = angle_ply(45, 4, 0.1, mat_cfrp)
    @test length(lam_angle) == 4
    @test isapprox(Laminates.thickness(lam_angle), 0.4, atol=1e-15)
end

@testset "Reference CLT and transformation checks" begin
    mat = OrthotropicMaterial(135e9, 10e9, 0.3, 5e9)
    q = Q12(mat)
    s = S12(mat)
    q_ref = q12_reference(mat)
    s_ref = s12_reference(mat)
    @test q ≈ q_ref rtol = 1e-12
    @test s ≈ s_ref rtol = 1e-12
    @test q * s ≈ identity_matrix(3) rtol = 1e-12 atol = 1e-12
    @test inv(q) ≈ s rtol = 1e-12 atol = 1e-12

    for θ in (0.0, 90.0, 45.0, -45.0)
        @test Qxy(Sheet(θ, 1.0, mat)) ≈ qbar_reference(mat, θ) rtol = 1e-10 atol = 1e-4
    end

    q45 = Qxy(Sheet(45.0, 1.0, mat))
    qm45 = Qxy(Sheet(-45.0, 1.0, mat))
    @test q45[1, 1] ≈ qm45[1, 1] rtol = 1e-12
    @test q45[2, 2] ≈ qm45[2, 2] rtol = 1e-12
    @test q45[1, 2] ≈ qm45[1, 2] rtol = 1e-12
    @test q45[3, 3] ≈ qm45[3, 3] rtol = 1e-12
    @test q45[1, 3] ≈ -qm45[1, 3] rtol = 1e-12 atol = 1e-6
    @test q45[2, 3] ≈ -qm45[2, 3] rtol = 1e-12 atol = 1e-6
end

@testset "Reference ABD integration and coupling behavior" begin
    mat = OrthotropicMaterial(130e9, 9e9, 0.28, 4.8e9)
    t = 0.125e-3

    lam_single = Laminate([0.0], [t], mat)
    q0 = qbar_reference(mat, 0.0)
    @test A(lam_single) ≈ q0 * t rtol = 1e-12
    @test B(lam_single) ≈ zeros(3, 3) atol = 1e-12
    @test D(lam_single) ≈ q0 * t^3 / 12 rtol = 1e-12

    θs = [0.0, 90.0]
    ts = [t, t]
    lam_2ply = Laminate(θs, ts, mat)
    a_ref, b_ref, d_ref = abd_reference(θs, ts, mat)
    @test A(lam_2ply) ≈ a_ref rtol = 1e-10 atol = 1e-4
    @test B(lam_2ply) ≈ b_ref rtol = 1e-10 atol = 1e-6
    @test D(lam_2ply) ≈ d_ref rtol = 1e-10 atol = 1e-10
    @test ABD(lam_2ply)[1:3, 1:3] ≈ a_ref rtol = 1e-10 atol = 1e-4
    @test ABD(lam_2ply)[1:3, 4:6] ≈ b_ref rtol = 1e-10 atol = 1e-6
    @test ABD(lam_2ply)[4:6, 4:6] ≈ d_ref rtol = 1e-10 atol = 1e-10

    lam_symmetric = Laminate([0.0, 90.0], t, mat; sym=:S)
    @test B(lam_symmetric) ≈ zeros(3, 3) atol = 1e-8

    lam_balanced = angle_ply(30.0, 4, t, mat)
    @test A(lam_balanced)[1, 3] ≈ 0.0 atol = 1e-6
    @test A(lam_balanced)[2, 3] ≈ 0.0 atol = 1e-6

    lam_cross = cross_ply(4, t, mat)
    @test D(lam_cross)[1, 3] ≈ 0.0 atol = 1e-8
    @test D(lam_cross)[2, 3] ≈ 0.0 atol = 1e-8
    @test !isapprox(D(lam_balanced)[1, 3], 0.0; atol=1e-8)
    @test !isapprox(D(lam_balanced)[2, 3], 0.0; atol=1e-8)

    lam_balanced_unsymmetric = Laminate([45.0, -45.0, 0.0], t, mat)
    @test !isapprox(D(lam_balanced_unsymmetric)[1, 3], 0.0; atol=1e-8)
    @test !isapprox(D(lam_balanced_unsymmetric)[2, 3], 0.0; atol=1e-8)
end

@testset "ABD inversion and load/deformation round-trip" begin
    mat = OrthotropicMaterial(140e9, 11e9, 0.27, 5.2e9)
    lam = Laminate([10.0, -35.0, 75.0], [0.10e-3, 0.12e-3, 0.08e-3], mat)
    abd_mtx = ABD(lam)
    abd_inv = abd(lam)
    @test abd_mtx * abd_inv ≈ identity_matrix(6) rtol = 1e-10 atol = 1e-10
    @test abd_inv * abd_mtx ≈ identity_matrix(6) rtol = 1e-10 atol = 1e-10

    load = [15_000.0, -9_000.0, 2_500.0, 1.8, -0.7, 0.9]
    ϵ0, κ = global_deformations(lam, load)
    deform = vcat(ϵ0, κ)
    @test abd_mtx * deform ≈ load rtol = 1e-10 atol = 1e-8
end

@testset "Ply-center and ply-surface stress recovery conventions" begin
    mat = OrthotropicMaterial(120e9, 8e9, 0.29, 4.3e9)
    θs = [30.0, -15.0]
    ts = [0.15e-3, 0.25e-3]
    lam = Laminate(θs, ts, mat)
    load = [8_000.0, -2_000.0, 1_500.0, 1.2, -0.6, 0.3]

    ϵ0, κ = global_deformations(lam, load)
    z = z_interfaces(ts)

    g_deforms_center = global_deformations_per_sheet(lam, load)
    g_tensions_center = global_tensions_per_sheet(lam, load)
    l_deforms_center = local_deformations(lam, load)
    l_tensions_center = local_tensions(lam, load)

    for i in eachindex(θs)
        zc = (z[i] + z[i + 1]) / 2
        ϵg_ref = ϵ0 .+ zc .* κ
        @test g_deforms_center[i] ≈ ϵg_ref rtol = 1e-12 atol = 1e-12

        qbar = qbar_reference(mat, θs[i])
        @test g_tensions_center[i] ≈ qbar * ϵg_ref rtol = 1e-9 atol = 1e-4
        @test l_deforms_center[i] ≈ Tϵ(θs[i]) * ϵg_ref rtol = 1e-12 atol = 1e-12
        @test l_tensions_center[i] ≈ Tσ(θs[i]) * g_tensions_center[i] rtol = 1e-10 atol =
            1e-4
    end

    for (i, θ) in pairs(θs), zsurf in (z[i], z[i + 1])
        ϵg = ϵ0 .+ zsurf .* κ
        σg = qbar_reference(mat, θ) * ϵg
        ϵl = Tϵ(θ) * ϵg
        σl = Tσ(θ) * σg
        @test σl ≈ q12_reference(mat) * ϵl rtol = 1e-10 atol = 1e-4
    end
end

@testset "Failure criteria reference checks" begin
    mat = OrthotropicMaterial(1.0, 1.0, 0.0, 1.0, 120.0, 90.0, 110.0, 80.0, 60.0)
    lam = Laminate([0.0], 1.0, mat)

    pure_load = [30.0, -20.0, 18.0, 0.0, 0.0, 0.0]
    @test max_tensions_criteria(lam, pure_load) == [[true, true, true]]
    @test max_strain_criteria(lam, pure_load) == [[true, true, true]]

    overload = [130.0, 95.0, 70.0, 0.0, 0.0, 0.0]
    @test max_tensions_criteria(lam, overload) == [[false, false, false]]
    @test max_strain_criteria(lam, overload) == [[false, false, false]]

    function tsai_hill_expected(σ1, σ2, σ6, m)
        if σ1 > 0 && σ2 > 0
            return (σ1 / m.Xt)^2 - (σ1 / m.Xt) * (σ2 / m.Xt) +
                   (σ2 / m.Yt)^2 +
                   (σ6 / m.S12)^2
        elseif σ1 < 0 && σ2 > 0
            return (σ1 / m.Xc)^2 +
                   (σ1 / m.Xc) * (σ2 / m.Xc) +
                   (σ2 / m.Yt)^2 +
                   (σ6 / m.S12)^2
        elseif σ1 > 0 && σ2 < 0
            return (σ1 / m.Xt)^2 +
                   (σ1 / m.Xt) * (σ2 / m.Xt) +
                   (σ2 / m.Yc)^2 +
                   (σ6 / m.S12)^2
        else
            return (σ1 / m.Xc)^2 - (σ1 / m.Xc) * (σ2 / m.Xc) +
                   (σ2 / m.Yc)^2 +
                   (σ6 / m.S12)^2
        end
    end

    for (σ1, σ2, σ6) in
        ((45.0, 25.0, 5.0), (-40.0, 25.0, -6.0), (40.0, -20.0, 7.0), (-35.0, -18.0, 4.0))
        fi = tsai_hill_criteria(lam, [σ1, σ2, σ6, 0.0, 0.0, 0.0])[1]
        @test fi ≈ tsai_hill_expected(σ1, σ2, σ6, mat) rtol = 1e-12 atol = 1e-12
    end

    function tsai_wu_expected(σ1, σ2, σ6, m)
        f1 = 1 / m.Xt - 1 / m.Xc
        f2 = 1 / m.Yt - 1 / m.Yc
        f11 = 1 / (m.Xt * m.Xc)
        f22 = 1 / (m.Yt * m.Yc)
        f66 = 1 / m.S12^2
        f12 = -sqrt(f11 * f22) / 2
        return (f1 * σ1 + f2 * σ2) +
               (f11 * σ1^2 + f22 * σ2^2 + f66 * σ6^2) +
               2 * f12 * σ1 * σ2
    end

    for (σ1, σ2, σ6) in ((30.0, 20.0, 10.0), (30.0, -20.0, 10.0), (-30.0, 20.0, -10.0))
        fi = tsai_wu_criteria(lam, [σ1, σ2, σ6, 0.0, 0.0, 0.0])[1]
        @test fi ≈ tsai_wu_expected(σ1, σ2, σ6, mat) rtol = 1e-12 atol = 1e-12
    end
end

@testset "Aqua" begin
    Aqua.test_all(Laminates; ambiguities=false, deps_compat=false)
end

@testset "Formatting" begin
    @test JuliaFormatter.format(Laminates; verbose=false, overwrite=false)
end
