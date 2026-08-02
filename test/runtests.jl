
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

# max_strain_criteria regression tests (issue #8)
# Material with E1=100, E2=10, G12=5, ν12=0.3.
# Strength allowables: Xt=500, Yt=50, Xc=300, Yc=30, S12=25.
# Derived strain allowables: ε1t=5.0, ε1c=3.0, ε2t=5.0, ε2c=3.0, γ12_allow=5.0.
# For a single 0° ply of thickness 1, applying load [Nx,Ny,Nxy,...]:
#   ε1 ≈ Nx/E1,  ε2 ≈ Ny/E2,  γ12 ≈ Nxy/G12  (Poisson coupling negligible for pure loads).
@testset "max_strain_criteria (issue #8)" begin
    mat = OrthotropicMaterial(100.0, 10.0, 0.3, 5.0, 500.0, 50.0, 300.0, 30.0, 25.0)
    lam = Laminates.Laminate([0.0], [1.0], mat)

    # Longitudinal tension: Nx=200 → ε1≈2.0 < ε1t=5.0 → pass
    @test Laminates.max_strain_criteria(lam, [200.0, 0, 0, 0, 0, 0])[1][1] == true
    # Longitudinal tension: Nx=600 → ε1≈6.0 > ε1t=5.0 → fail
    @test Laminates.max_strain_criteria(lam, [600.0, 0, 0, 0, 0, 0])[1][1] == false

    # Longitudinal compression: Nx=-100 → ε1≈-1.0 > -ε1c=-3.0 → pass
    @test Laminates.max_strain_criteria(lam, [-100.0, 0, 0, 0, 0, 0])[1][1] == true
    # Longitudinal compression: Nx=-400 → ε1≈-4.0 < -ε1c=-3.0 → fail
    @test Laminates.max_strain_criteria(lam, [-400.0, 0, 0, 0, 0, 0])[1][1] == false

    # Transverse tension: Ny=30 → ε2≈3.0 < ε2t=5.0 → pass
    @test Laminates.max_strain_criteria(lam, [0, 30.0, 0, 0, 0, 0])[1][2] == true
    # Transverse tension: Ny=60 → ε2≈6.0 > ε2t=5.0 → fail
    @test Laminates.max_strain_criteria(lam, [0, 60.0, 0, 0, 0, 0])[1][2] == false

    # Transverse compression: Ny=-20 → ε2≈-2.0 > -ε2c=-3.0 → pass
    @test Laminates.max_strain_criteria(lam, [0, -20.0, 0, 0, 0, 0])[1][2] == true
    # Transverse compression: Ny=-40 → ε2≈-4.0 < -ε2c=-3.0 → fail
    @test Laminates.max_strain_criteria(lam, [0, -40.0, 0, 0, 0, 0])[1][2] == false

    # Pure shear: Nxy=15 → γ12≈3.0 < γ12_allow=5.0 → pass
    @test Laminates.max_strain_criteria(lam, [0, 0, 15.0, 0, 0, 0])[1][3] == true
    # Pure shear: Nxy=30 → γ12≈6.0 > γ12_allow=5.0 → fail
    @test Laminates.max_strain_criteria(lam, [0, 0, 30.0, 0, 0, 0])[1][3] == false

    # Dimensional-correctness proof: with Xt=300 → ε1t=3.0.
    # Nx=400 → ε1≈4.0, which satisfies -300 < 4.0 < 300 (wrong stress comparison)
    # but fails the correct strain comparison -3.0 < 4.0 < 3.0.
    mat2 = OrthotropicMaterial(100.0, 10.0, 0.3, 5.0, 300.0, 50.0, 300.0, 30.0, 25.0)
    lam2 = Laminates.Laminate([0.0], [1.0], mat2)
    @test Laminates.max_strain_criteria(lam2, [400.0, 0, 0, 0, 0, 0])[1][1] == false

    # Combined load: Nx=200, Ny=20 both within allowables → all pass
    @test all(Laminates.max_strain_criteria(lam, [200.0, 20.0, 0, 0, 0, 0])[1])
    # Combined load: Nx=200, Ny=60 → transverse direction fails
    crit_combined = Laminates.max_strain_criteria(lam, [200.0, 60.0, 0, 0, 0, 0])[1]
    @test crit_combined[1] == true
    @test crit_combined[2] == false
end

@testset "Aqua" begin
    Aqua.test_all(Laminates; ambiguities=false, deps_compat=false)
end

@testset "Formatting" begin
    @test JuliaFormatter.format(Laminates; verbose=false, overwrite=false)
end
