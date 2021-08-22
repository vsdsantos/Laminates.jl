
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
