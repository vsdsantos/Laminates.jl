
using LinearAlgebra

struct OrthotropicMaterial{T<:Real}
    E1::T
    E2::T
    ν12::T
    G12::T
    Xt::Union{T,Nothing}
    Yt::Union{T,Nothing}
    Xc::Union{T,Nothing}
    Yc::Union{T,Nothing}
    S12::Union{T,Nothing}

    function OrthotropicMaterial(E1::T, E2::T, nu12::T, G12::T) where {T<:Real}
        return new{T}(E1, E2, nu12, G12, nothing, nothing, nothing, nothing, nothing)
    end

    function OrthotropicMaterial(
        E1::T, E2::T, nu12::T, G12::T, Xt::T, Yt::T, Xc::T, Yc::T, S12::T
    ) where {T<:Real}
        return new{T}(E1, E2, nu12, G12, Xt, Yt, Xc, Yc, S12)
    end
end

function OrthotropicMaterial(E1, E2, nu12, G12)
    T = promote_type(typeof(E1), typeof(E2), typeof(nu12), typeof(G12))
    return OrthotropicMaterial(T(E1), T(E2), T(nu12), T(G12))
end

function OrthotropicMaterial(E1, E2, nu12, G12, Xt, Yt, Xc, Yc, S12)
    T = promote_type(
        typeof(E1),
        typeof(E2),
        typeof(nu12),
        typeof(G12),
        typeof(Xt),
        typeof(Yt),
        typeof(Xc),
        typeof(Yc),
        typeof(S12),
    )
    return OrthotropicMaterial(
        T(E1), T(E2), T(nu12), T(G12), T(Xt), T(Yt), T(Xc), T(Yc), T(S12)
    )
end

function Base.show(io::IO, mat::OrthotropicMaterial)
    return print(io, "E1=$(mat.E1),E2=$(mat.E2),ν12=$(mat.ν12),G12=$(mat.G12)")
end

function Q12(mat::OrthotropicMaterial)
    Q11 = mat.E1^2/(mat.E1-mat.ν12^2*mat.E2)
    Q12 = mat.ν12*mat.E1*mat.E2/(mat.E1-mat.ν12^2*mat.E2)
    Q22 = mat.E1*mat.E2/(mat.E1-mat.ν12^2*mat.E2)
    Q66 = mat.G12
    return [
        Q11 Q12 0
        Q12 Q22 0
        0 0 Q66
    ]
end

function S12(mat::OrthotropicMaterial)
    S11 = 1/mat.E1
    S12 = -mat.ν12/mat.E1
    S22 = 1/mat.E2
    S66 = 1/mat.G12
    return [
        S11 S12 0
        S12 S22 0
        0 0 S66
    ]
end
