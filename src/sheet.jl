
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