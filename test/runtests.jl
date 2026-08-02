
using Laminates
using Laminates: OrthotropicMaterial, cross_ply, angle_ply
using Test
using Aqua
using JuliaFormatter

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

@testset "Tsai-Hill criteria" begin
    mat = OrthotropicMaterial(130e9, 9e9, 0.28, 5e9, 1000.0, 40.0, 800.0, 200.0, 80.0)
    lam = Laminates.Laminate([0.0], 1.0, mat)

    load_from_local_stress(σ1, σ2, τ12) = [σ1, σ2, τ12, 0.0, 0.0, 0.0]

    # Pure local stress components
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(500.0, 0.0, 0.0))[1],
        (500 / 1000)^2;
        atol=1e-10,
    )
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(0.0, 20.0, 0.0))[1],
        (20 / 40)^2;
        atol=1e-10,
    )
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(0.0, 0.0, 40.0))[1],
        (40 / 80)^2;
        atol=1e-10,
    )

    # All σ1/σ2 sign quadrants
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(500.0, 20.0, 0.0))[1],
        (500 / 1000)^2 - (500 * 20) / (1000 * 40) + (20 / 40)^2;
        atol=1e-10,
    )
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(-400.0, 20.0, 0.0))[1],
        (-400 / 800)^2 - (-400 * 20) / (800 * 40) + (20 / 40)^2;
        atol=1e-10,
    )
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(500.0, -60.0, 0.0))[1],
        (500 / 1000)^2 - (500 * -60) / (1000 * 200) + (-60 / 200)^2;
        atol=1e-10,
    )
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(-400.0, -60.0, 0.0))[1],
        (-400 / 800)^2 - (-400 * -60) / (800 * 200) + (-60 / 200)^2;
        atol=1e-10,
    )

    # Strongly orthotropic regression case for the interaction denominator.
    # Correct: -σ1σ2/(Xt*Yt) = -(300*8)/(1000*40) = -0.06 (old Xt*Xt bug gives -0.0024).
    @test isapprox(
        Laminates.tsai_hill_criteria(lam, load_from_local_stress(300.0, 8.0, 0.0))[1],
        0.09 - 0.06 + 0.04;
        atol=1e-10,
    )
end

@testset "Aqua" begin
    Aqua.test_all(Laminates; ambiguities=false, deps_compat=false)
end

@testset "Formatting" begin
    @test JuliaFormatter.format(Laminates; verbose=false, overwrite=false)
end
