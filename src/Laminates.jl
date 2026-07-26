module Laminates

include("materials.jl")
include("sheet.jl")

struct Laminate
    sheets::AbstractVector{<:Sheet}
end

function symmetry_operations(vec::AbstractVector, sym::Symbol)
    if sym == :S
        return vcat(vec, reverse(deepcopy(vec)))
    elseif sym == :SM
        return vcat(vec, reverse(deepcopy(vec[1:(end - 1)])))
    elseif sym == :AS
        return vcat(vec, deepcopy(vec))
    elseif sym == :ASM
        return vcat(vec[1:end], deepcopy(vec[1:(end - 1)]))
    else
        return vec
    end
end

function Laminate(
    θs::AbstractVector{<:Real},
    ts::AbstractVector{<:Real},
    mats::AbstractVector{<:OrthotropicMaterial};
    sym=:no,
)
    N = length(θs)
    if length(ts) != N || length(mats) != N
        error("Length of arrays are not the same.")
    end
    return Laminate(symmetry_operations([Sheet(θs[i], ts[i], mats[i]) for i in 1:N], sym))
end

function Laminate(
    θs::AbstractVector{<:Real},
    t::Real,
    mats::AbstractVector{<:OrthotropicMaterial};
    sym=:no,
)
    N = length(θs)
    if length(mats) != N
        error("Length of arrays are not the same.")
    end
    return Laminate(symmetry_operations([Sheet(θs[i], t, mats[i]) for i in 1:N], sym))
end

function Laminate(
    θs::AbstractVector{<:Real},
    ts::AbstractVector{<:Real},
    mat::OrthotropicMaterial;
    sym=:no,
)
    N = length(θs)
    if length(ts) != N
        error("Length of arrays are not the same.")
    end
    return Laminate(symmetry_operations([Sheet(θs[i], ts[i], mat) for i in 1:N], sym))
end

function Laminate(θs::AbstractVector{<:Real}, t::Real, mat::OrthotropicMaterial; sym=:no)
    N = length(θs)
    return Laminate(symmetry_operations([Sheet(θs[i], t, mat) for i in 1:N], sym))
end

function angle_ply(θ::Real, N::Integer, t::Real, mat::OrthotropicMaterial)::Laminate
    if N % 2 != 0
        error("N should be even.")
    end
    thetas = [isodd(i) ? θ : -θ for i in 1:(N ÷ 2)]
    return Laminate(thetas, t, mat; sym=:S)
end

function cross_ply(N::Integer, t::Real, mat::OrthotropicMaterial)::Laminate
    thetas = [iseven(i) ? 90 : 0 for i in 1:(Int64(ceil(N / 2)))]
    if N % 2 != 0
        return Laminate(thetas, t, mat; sym=:SM)
    else
        return Laminate(thetas, t, mat; sym=:S)
    end
end

thickness(lam::Laminate) = sum(sh.t_i for sh in lam.sheets)

Base.length(lam::Laminate) = Base.length(lam.sheets)

function Base.show(io::IO, lam::Laminate)
    results = [(l.θ_i, l.t_i) for l in lam.sheets]
    return print(io, "(θ,t)=$(results)")
end

function t_pos(lam::Laminate)
    t = thickness(lam)
    t_k = -t / 2
    ts = Float64[]
    for sh in lam.sheets
        push!(ts, t_k)
        t_k += sh.t_i
    end
    push!(ts, t / 2)
    return ts
end

function A(lam::Laminate)
    t_k = t_pos(lam)
    result = zeros(3, 3)
    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        result += Qxy(sh) * (t_k[i + 1] - t_k[i])
    end
    return result
end

function B(lam::Laminate)
    t_k = t_pos(lam)
    result = zeros(3, 3)
    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        result += Qxy(sh) * (t_k[i + 1]^2 - t_k[i]^2)
    end
    return result ./ 2
end

function D(lam::Laminate)
    t_k = t_pos(lam)
    result = zeros(3, 3)
    for i in eachindex(lam.sheets)
        sh = lam.sheets[i]
        result += Qxy(sh) * (t_k[i + 1]^3 - t_k[i]^3)
    end
    return result ./ 3
end

function ABD(lam::Laminate)
    A_, B_, D_ = A(lam), B(lam), D(lam)
    return [
        A_ B_
        B_ D_
    ]
end

function abd(lam::Laminate)
    return inv(ABD(lam))
end

function equivalent_material(lam::Laminate)
    abd_ = abd(lam)
    a = abd_[1:3, 1:3]
    t = thickness(lam)

    Ex = 1 / (t * a[1, 1])
    Ey = 1 / (t * a[2, 2])
    Gxy = 1 / (t * a[3, 3])
    νxy = a[2, 1] / a[1, 1]

    return OrthotropicMaterial(Ex, Ey, νxy, Gxy)
end

include("strain_stress.jl")
include("failure_criterias.jl")

function Tσ(θ::Real)
    θ = θ * π / 180
    m, n = cos(θ), sin(θ)

    return [
        m^2 n^2 2*m*n
        n^2 m^2 -2*m*n
        -m*n m*n m^2-n^2
    ]
end

function Tϵ(θ::Real)
    θ = θ * π / 180
    m, n = cos(θ), sin(θ)

    return [
        m^2 n^2 m*n
        n^2 m^2 -m*n
        -2*m*n 2*m*n m^2-n^2
    ]
end

end # module
