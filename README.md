# Laminates.jl
Tools for modeling composite laminates.

## Using

```julia
using Laminates: OrthotropicMaterial, cross_ply, ABD, tsai_wu_criteria

cfrp = OrthotropicMaterial(70e9, 20e9, 0.3, 10e8, 100,100,100,100,20)

lam = cross_ply(5, 0.1, cfrp)

ABD(lam)
```
```
6×6 Matrix{Float64}:
  2.56598e10   3.07918e9    -1.51483e-7   -2.38419e-7   0.0          -1.65436e-24
  3.07918e9    2.05279e10    7.79967e-7    0.0          0.0           6.61744e-24
 -1.51483e-7   7.79967e-7    5.0e8        -1.65436e-24  6.61744e-24  -3.72529e-9
 -2.38419e-7   0.0          -1.65436e-24   6.37219e8    6.41496e7    -1.64106e-9
  0.0          0.0           6.61744e-24   6.41496e7    3.25024e8     8.44964e-9
 -1.65436e-24  6.61744e-24  -3.72529e-9   -1.64106e-9   8.44964e-9    1.04167e7
```
```julia
tsai_wu_criteria(lam, [100,100,10,0,0,0])
```
```
5-element Vector{Any}:
 6.495762834897815
 9.892961286870392
 6.495762834897815
 9.892961286870392
 6.495762834897818
```
