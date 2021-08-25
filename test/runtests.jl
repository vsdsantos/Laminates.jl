
using Laminates
using Laminates: OrthotropicMaterial, cross_ply, angle_ply
using Test

mat_cfrp = OrthotropicMaterial(70e9, 20e9, 0.3, 10e8)
mat_cfrp_limit = OrthotropicMaterial(70e9, 20e9, 0.3, 10e8, 100,100,100,100,20)

lam_cross = cross_ply(10, 0.1, mat_cfrp)

@test length(lam_cross) == 10

@test isapprox(Laminates.thickness(lam_cross), 1.0, atol=1e-15)