
struct Sheet{T<:Real,M<:OrthotropicMaterial}
    θ_i::T
    t_i::T
    material::M

    function Sheet(θ_i::T, t_i::T, material::M) where {T<:Real,M<:OrthotropicMaterial}
        return new{T,M}(θ_i, t_i, material)
    end
end

function Sheet(θ_i, t_i, material::OrthotropicMaterial)
    T = promote_type(typeof(θ_i), typeof(t_i))
    return Sheet(T(θ_i), T(t_i), material)
end

function Qxy(sh::Sheet)
    Q_local = Q12(sh.material)
    Tsig = Tσ(-sh.θ_i) # matriz de transformação inversa
    Teps = Tϵ(sh.θ_i) # matriz de transformação inversa
    return Tsig * Q_local * Teps
end
