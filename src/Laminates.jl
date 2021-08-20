module Laminates

using LinearAlgebra

struct OrthotropicMaterial
	E1::Real
	E2::Real
	ν12::Real
	G12::Real
	Xt::Real
	Yt::Real
	Xc::Real
	Yc::Real
	S12::Real
	
	OrthotropicMaterial(E1, E2, nu12, G12) = new(E1, E2, nu12, G12)
	OrthotropicMaterial(E1, E2, nu12, G12, Xt, Yt, Xc, Yc, S12) = new(E1, E2, nu12, G12, Xt, Yt, Xc, Yc, S12)
end

function Q12(mat::OrthotropicMaterial)
	Q11 = mat.E1^2/(mat.E1-mat.ν12^2*mat.E2)
	Q12 = mat.ν12*mat.E1*mat.E2/(mat.E1-mat.ν12^2*mat.E2)
	Q22 = mat.E1*mat.E2/(mat.E1-mat.ν12^2*mat.E2)
	Q66 = mat.G12
	[Q11 Q12 0
	 Q12 Q22 0
	  0   0  Q66]
end

function S12(mat::OrthotropicMaterial)
	S11 = 1/mat.E1
	S12 = -mat.ν12/mat.E1
	S22 = 1/mat.E2
	S66 = 1/mat.G12
	[S11 S12  0
	 S12 S22  0
	 0   0  S66]
end

struct Sheet
	θ_i::Real
	t_i::Real
	material::OrthotropicMaterial
end

function Qxy(sh::Sheet)
	Q_local = Q12(sh.material)
	Tsig = Tσ(-sh.θ_i) # matriz de transformação inversa
	Teps = Tϵ(sh.θ_i) # matriz de transformação inversa
	Tsig * Q_local * Teps
end

struct Laminate
	sheets::Array{Sheet}
end

thickness(lam::Laminate) = sum([sh.t_i for sh in lam.sheets])

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
	for i in 1:length(lam.sheets)
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

function global_deformations(lam::Laminate, load)
	deform = ABD(lam) \ load
	deform[1:3], deform[4:6]
end

function global_deformations_per_sheet(lam::Laminate, load)
	ϵ, κ = global_deformations(lam, load)
	t_k = t_pos(lam)
	deforms = []
	for i in 1:length(lam.sheets)
		sh = lam.sheets[i]
		push!(deforms, ϵ .+ (t_k[i]+sh.t_i/2) .*κ)
	end
	deforms
end

function global_tensions_per_sheet(lam::Laminate, load)
	g_deforms = global_deformations_per_sheet(lam, load)
	
	tensions = []
	
	for i in 1:length(lam.sheets)
		sh = lam.sheets[i]
		push!(tensions, Qxy(sh) * g_deforms[i])
	end
	
	tensions
end

function local_deformations(lam::Laminate, load)
	g_deforms = global_deformations_per_sheet(lam, load)
	
	deforms = []
	
	for i in 1:length(lam.sheets)
		sh = lam.sheets[i]
		
		push!(deforms, Tϵ(sh.θ_i) * g_deforms[i])
	end
	
	deforms
end

function local_tensions(lam::Laminate, load)
	g_tensions = global_tensions_per_sheet(lam, load)
	
	tensions = []
	
	for i in 1:length(lam.sheets)
		sh = lam.sheets[i]
		push!(tensions, Tσ(sh.θ_i) * g_tensions[i])
	end
	
	tensions
end

function max_tensions_criteria(lam::Laminate, load)
	tensions = local_tensions(lam, load)
	
	crit = []
	
	for i in 1:length(lam.sheets)
		m = lam.sheets[i].material
		σ = tensions[i]
		push!(crit, [-m.Xc < σ[1] < m.Xt, -m.Yc < σ[2] < m.Yt, abs(σ[3]) < m.S12])
	end
	
	crit
end

function max_strain_criteria(lam::Laminate, load)
	strain = local_deformations(lam, load)
	
	crit = []
	
	for i in 1:length(lam.sheets)
		m = lam.sheets[i].material
		ϵ = strain[i]
		push!(crit, [-m.Xc < ϵ[1] < m.Xt, -m.Yc < ϵ[2] < m.Yt, abs(ϵ[3]) < m.S12])
	end
	
	crit
end

function tsai_hill_criteria(lam::Laminate, load)
	tensions = local_tensions(lam, load)
	
	crit = []
	
	for i in 1:length(lam.sheets)
		m = lam.sheets[i].material
		σ = tensions[i]
		σ1, σ2, σ3 = σ
		if σ1 > 0 && σ2 > 0
			FI = (σ1/m.Xt)^2 - (σ1/m.Xt)*(σ2/m.Xt) + (σ2/m.Yt)^2 + (σ3/m.S12)^2
		elseif σ1 < 0 && σ2 > 0
			FI = (σ1/m.Xc)^2 + (σ1/m.Xc)*(σ2/m.Xc) + (σ2/m.Yt)^2 + (σ3/m.S12)^2
		elseif σ1 > 0 && σ2 < 0
			FI = (σ1/m.Xt)^2 + (σ1/m.Xt)*(σ2/m.Xt) + (σ2/m.Yc)^2 + (σ3/m.S12)^2
		else σ1 < 0 && σ2 < 0
			FI = (σ1/m.Xc)^2 - (σ1/m.Xc)*(σ2/m.Xc) + (σ2/m.Yc)^2 + (σ3/m.S12)^2
		end
		push!(crit, FI)
	end
	
	crit
end

function tsai_wu_criteria(lam::Laminate, load)
	tensions = local_tensions(lam, load)
	
	crit = []
	
	for i in 1:length(lam.sheets)
		m = lam.sheets[i].material
		
		F1 = 1/m.Xt - 1/m.Xc
		F2 = 1/m.Yt - 1/m.Yc
		F11 = 1 / (m.Xt * m.Xc)
		F22 = 1 / (m.Yt * m.Yc)
		F66 = 1 / (m.S12^2)
		F12 = -sqrt(F11*F22)/2
		
		σ = tensions[i]
		σ1, σ2, σ6 = σ
		
		FI = (F1*σ1+F2*σ2)+(F11*σ1^2+F22*σ2^2+F66*σ6^2)+2*F12*σ1*σ2
		
		push!(crit, FI)
	end
	
	crit
end

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
