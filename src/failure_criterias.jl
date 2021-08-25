
function max_tensions_criteria(lam::Laminate, load::AbstractVector{<:Real})
	tensions = local_tensions(lam, load)
	
	crit = []
	
	for i in 1:length(lam.sheets)
		m = lam.sheets[i].material
		σ = tensions[i]
		push!(crit, [-m.Xc < σ[1] < m.Xt, -m.Yc < σ[2] < m.Yt, abs(σ[3]) < m.S12])
	end
	
	crit
end

function max_strain_criteria(lam::Laminate, load::AbstractVector{<:Real})
	strain = local_deformations(lam, load)
	
	crit = []
	
	for i in 1:length(lam)
		m = lam.sheets[i].material
		ϵ = strain[i]
		push!(crit, [-m.Xc < ϵ[1] < m.Xt, -m.Yc < ϵ[2] < m.Yt, abs(ϵ[3]) < m.S12])
	end
	
	crit
end

function tsai_hill_criteria(lam::Laminate, load::AbstractVector{<:Real})
	tensions = local_tensions(lam, load)
	
	crit = []
	
	for i in 1:length(lam)
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

function tsai_wu_criteria(lam::Laminate, load::AbstractVector{<:Real})
	tensions = local_tensions(lam, load)
	
	crit = []
	
	for i in 1:length(lam)
		m = lam.sheets[i].material
		
		if m.Xt == undef || m.Yt == undef || m.Xc == undef || m.Yc == undef || m.S12 == undef
			error("Some Material properties are undefined.")
		end

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