module Laminates

import Base: length

include("materials.jl")
include("sheet.jl")

struct Laminate
	sheets::AbstractVector{<:Sheet}
end

function symetry_operations(vec::AbstractVector, sym::Symbol)
	if sym == :S
		return vcat(vec, reverse(vec))
	elseif sym == :SM
		return vcat(vec[1:end-1], [vec[end]], reverse(vec[1:end-1]))
	elseif sym == :AS
		return vcat(vec, vec)
	elseif sym == :ASM
		return vcat(vec[1:end-1], [vec[end]], vec[1:end-1])
	else
		return vec
	end
end

function Laminate(θs::AbstractVector{<:Real}, ts::AbstractVector{<:Real}, mats::AbstractVector{<:OrthotropicMaterial}; sym = :no)
	N = length(θs)
	if length(ts) != N || length(mats) != N
		error("Length of arrays are not the same.")
	end
	Laminate(symetry_operations([Sheet(θs[i], ts[i], mats[i]) for i in 1:N], sym))
end

function Laminate(θs::AbstractVector{<:Real}, t::Real, mats::AbstractVector{<:OrthotropicMaterial}; sym = :no)
	N = length(θs)
	if length(mats) != N
		error("Length of arrays are not the same.")
	end
	Laminate(symetry_operations([Sheet(θs[i], t, mats[i]) for i in 1:N], sym))
end

function Laminate(θs::AbstractVector{<:Real}, ts::AbstractVector{<:Real}, mat::OrthotropicMaterial; sym = :no)
	N = length(θs)
	if length(ts) != N
		error("Length of arrays are not the same.")
	end
	Laminate(symetry_operations([Sheet(θs[i], ts[i], mat) for i in 1:N], sym))
end

function Laminate(θs::AbstractVector{<:Real}, t::Real, mat::OrthotropicMaterial; sym = :no)
	N = length(θs)
	Laminate(symetry_operations([Sheet(θs[i], t, mat) for i in 1:N], sym))
end

function angle_ply(θ::Real, N::Integer, t::Real,  mat::OrthotropicMaterial)::Laminate
	if N % 2 != 0
		error("N should be even.")
	end
	thetas = Vector{Real}()
	for i in 1:(N÷2)
		if i % 2 == 0
			push!(thetas, -θ)
		else
			push!(thetas, θ)
		end
	end
	Laminate(thetas, t, mat, sym=:S)
end

function cross_ply(N::Integer, t::Real,  mat::OrthotropicMaterial)::Laminate
	if N % 2 != 0
		sym = :SM
	else
		sym = :S
	end
	thetas = Vector{Real}()
	for i in 1:(N÷2)
		if i % 2 == 0
			push!(thetas, 90)
		else
			push!(thetas, 0)
		end
	end
	Laminate(thetas, t, mat, sym=sym)
end

thickness(lam::Laminate) = sum([sh.t_i for sh in lam.sheets])

length(lam::Laminate) = length(lam.sheets)

function t_pos(lam::Laminate)
	t = thickness(lam)
	t_k = -t/2
	ts = []
	for sh in lam.sheets
		push!(ts, t_k)
		t_k += sh.t_i
	end
	push!(ts, t/2)
	ts
end

function A(lam::Laminate)
	t_k = t_pos(lam)
	A = zeros(3,3)
	for i in 1:length(lam)
		sh = lam.sheets[i]
		A += Qxy(sh) * (t_k[i+1]-t_k[i])
	end
	A
end

function B(lam::Laminate)
	t_k = t_pos(lam)
	B = zeros(3,3)
	for i in 1:length(lam.sheets)
		sh = lam.sheets[i]
		B += Qxy(sh) * (t_k[i+1]^2-t_k[i]^2)
	end
	B ./ 2
end

function D(lam::Laminate)
	t_k = t_pos(lam)
	D = zeros(3,3)
	for i in 1:length(lam.sheets)
		sh = lam.sheets[i]
		D += Qxy(sh) * (t_k[i+1]^3-t_k[i]^3)
	end
	D ./ 3
end

function ABD(lam::Laminate)
	A_, B_, D_ = A(lam), B(lam), D(lam)
	[A_ B_
	 B_ D_]
end

function abd(lam::Laminate)
	inv(ABD(lam))
end

function equivalent_material(lam::Laminate)
	abd_ = abd(lam)
	a = abd_[1:3,1:3]
	t = thickness(lam)
	
	Ex = 1/(t*a[1,1])
	Ey = 1/(t*a[2,2])
	Gxy = 1/(t*a[3,3])
	νxy = a[2,1]/a[1,1]
	
	OrthotropicMaterial(Ex, Ey, νxy, Gxy)
end

include("strain_stress.jl")
include("failure_criterias.jl")

function Tσ(θ::Real)
	θ = θ * π / 180
	m, n = cos(θ), sin(θ)
	
	[m^2  n^2   2*m*n
	 n^2  m^2  -2*m*n
	-m*n  m*n  m^2-n^2]
end

function Tϵ(θ::Real)
	θ = θ * π / 180
	m, n = cos(θ), sin(θ)
	
	[ m^2   n^2    m*n
	  n^2   m^2   -m*n
	-2*m*n 2*m*n m^2-n^2]
end

end # module
