import ForwardDiff: Dual

# TODO Based on this type do some generic things like use spin-scaling relations
# Note: Kind is needed because GGA enhancement or spin scaling etc. work differently
#       for exchange and correlation ... if it does not help, remove it again
abstract type Functional{Family,Kind} end

"""Return the family of a functional. Results are `:lda`, `:gga`, `:mgga` and
`:mggal` (Meta-GGA requiring Laplacian of ρ)"""
family(::Functional{F}) where {F} = F

"""
Return the functional kind: `:x` (exchange), `:c` (correlation), `:k` (kinetic) or
`:xc` (exchange and correlation combined)
"""
kind(::Functional{F,K}) where {F,K} = K

"""
Return the identifier corresponding to a functional, if available,
and `nothing` otherwise.
"""
function identifier end
function Base.show(io::IO, fun::Functional)
    id = identifier(fun)
    if isnothing(id)
        Base.show_default(io, fun)
    else
        print(io, id)
    end
end

@doc raw"""
True if the functional needs ``σ = 𝛁ρ ⋅ 𝛁ρ``.
"""
needs_σ(::Functional{F}) where {F} = (F in (:gga, :mgga, :mggal))

@doc raw"""
True if the functional needs ``τ`` (kinetic energy density).
"""
needs_τ(::Functional{F}) where {F} = (F in (:mgga, :mggal))

@doc raw"""
True if the functional needs ``Δ ρ``.
"""
needs_Δρ(::Functional{F}) where {F} = (F in (:mggal,))

"""
Does this functional support energy evaluations? Some don't, in which case
energy terms will not be returned by `potential_terms` and `kernel_terms`,
i.e. `e` will be `false` (a strong zero).
"""
has_energy(::Functional) = true

"""
The type of the parameters of the functional. If the functional has multiple parameters,
the result of `promote_type(...)` of the parameter types should be returned.

This allows the functional to influence the computed type, for example if it contains
Dual numbers.
"""
parameters_type(::Functional) = Bool # default: Bool, which promotes to any number type

# TODO These values are read-only for now and their defaults hard-coded for Float64
"""
Threshold for the density (below this value, functionals and derivatives
evaluate to zero). The threshold may depend on the floating-point type used
to represent densities and potentials, which is passed as the second argument.
"""
threshold_ρ(::Functional, T=Float64) = T(1e-15)  # TODO This might differ between functionals
threshold_σ(f::Functional, T=Float64) = threshold_ρ(f, T)^(4 // 3)
threshold_τ(::Functional, T=Float64)  = T(1e-20)
threshold_ζ(::Functional, T=Float64)  = eps(T)

# Drop dual types from threshold functions
threshold_ρ(f::Functional, T::Type{<:Dual}) = threshold_ρ(f, ForwardDiff.valtype(T))
threshold_σ(f::Functional, T::Type{<:Dual}) = threshold_σ(f, ForwardDiff.valtype(T))
threshold_τ(f::Functional, T::Type{<:Dual}) = threshold_τ(f, ForwardDiff.valtype(T))
threshold_ζ(f::Functional, T::Type{<:Dual}) = threshold_ζ(f, ForwardDiff.valtype(T))

# Silently drop extra arguments from evaluation functions
for fun in (:potential_terms, :kernel_terms)
    @eval begin
        $fun(func::Functional{:lda}, ρ, σ, args...)         = $fun(func, ρ)
        $fun(func::Functional{:gga}, ρ, σ, τ, args...)      = $fun(func, ρ, σ)
        $fun(func::Functional{:mgga}, ρ, σ, τ, Δρ, args...) = $fun(func, ρ, σ, τ)
    end
end

@doc raw"""
    potential_terms(f::Functional, ρ, [σ, τ, Δρ])

Evaluate energy and potential terms at a real-space grid of densities, density
derivatives etc. Not required derivatives for the functional type will be ignored.
Returns a named tuple with keys `e` (Energy per unit volume),
`Vρ` (``\frac{∂e}{∂ρ}``), `Vσ` (``\frac{∂e}{∂σ}``),
`Vτ` (``\frac{∂e}{∂τ}``), `Vl` (``\frac{∂e}{∂(Δρ)}``).
"""
function potential_terms end

@doc raw"""
    kernel_terms(f::Functional, ρ, [σ, τ, Δρ])

Evaluate energy, potential and kernel terms at a real-space grid of densities, density
derivatives etc. Not required derivatives for the functional type will be ignored.
Returns a named tuple with the same keys as `potential_terms` and additionally
second-derivative cross terms such as `Vρσ` (``\frac{∂^2e}{∂ρ∂σ}``).
"""
function kernel_terms end

#
# LDA
#
function potential_terms(func::Functional{:lda}, ρ::AbstractMatrix{T}) where {T}
    @assert has_energy(func)  # Otherwise custom implementation of this function needed
    s_ρ, n_p = size(ρ)
    TT = arithmetic_type(func, T)

    e  = similar(ρ, TT, n_p)
    Vρ = similar(ρ, TT, s_ρ, n_p)
    @views for i = 1:n_p
        potential_terms!(e[i:i], Vρ[:, i], func, ρ[:, i])
    end
    (; e, Vρ)
end
function potential_terms!(e, Vρ, func::Functional{:lda}, ρ::AbstractVector{T}) where {T}
    res = ForwardDiff.gradient!(DiffResults.DiffResult(zero(eltype(e)), Vρ),
                                ρ -> energy(func, ρ), ρ)
    e .= DiffResults.value(res)
    nothing
end

function kernel_terms(func::Functional{:lda}, ρ::AbstractMatrix{T}) where {T}
    @assert has_energy(func)
    s_ρ, n_p = size(ρ)
    TT = arithmetic_type(func, T)

    e   = similar(ρ, TT, n_p)
    Vρ  = similar(ρ, TT, s_ρ, n_p)
    Vρρ = similar(ρ, TT, s_ρ, s_ρ, n_p)

    # TODO Needed to make forward-diff work with !isbits floating-point types (e.g. BigFloat)
    Vρ  .= zero(T)
    Vρρ .= zero(T)

    @views for i = 1:n_p
        kernel_terms!(e[i:i], Vρ[:, i], Vρρ[:, :, i], func, ρ[:, i])
    end
    (; e, Vρ, Vρρ)
end
function kernel_terms!(e, Vρ, Vρρ, func::Functional{:lda}, ρ::AbstractVector{T}) where {T}
    res = ForwardDiff.hessian!(DiffResults.DiffResult(zero(eltype(e)), Vρ, Vρρ),
                               ρ -> energy(func, ρ), ρ)
    e .= DiffResults.value(res)
    nothing
end

function energy(func::Functional{:lda}, ρ::AbstractVector{T}) where {T}
    length(ρ) == 1 || error("Multiple spins not yet implemented for fallback functionals")
    ρtotal = ρ[1]
    if ρtotal ≤ threshold_ρ(func, T)
        zero(T)
    else
        energy(func, ρtotal)
    end
end

#
# GGA
#
function potential_terms(func::Functional{:gga}, ρ::AbstractMatrix{T},
                         σ::AbstractMatrix{U}) where {T,U}
    @assert has_energy(func)  # Otherwise custom implementation of this function needed
    s_ρ, n_p = size(ρ)
    s_σ = size(σ, 1)
    TT = arithmetic_type(func, T, U)

    e  = similar(ρ, TT, n_p)
    Vρ = similar(ρ, TT, s_ρ, n_p)
    Vσ = similar(ρ, TT, s_σ, n_p)
    @views for i = 1:n_p
        potential_terms!(e[i:i], Vρ[:, i], Vσ[:, i], func, ρ[:, i], σ[:, i])
    end
    (; e, Vρ, Vσ)
end
function potential_terms!(e, Vρ, Vσ, func::Functional{:gga},
                          ρ::AbstractVector, σ::AbstractVector)
    res = ForwardDiff.gradient!(DiffResults.DiffResult(zero(eltype(e)), Vρ),
                                ρ -> energy(func, ρ, σ), ρ)
    ForwardDiff.gradient!(DiffResults.DiffResult(zero(eltype(e)), Vσ),
                          σ -> energy(func, ρ, σ), σ)
    e .= DiffResults.value(res)
    nothing
end

function kernel_terms(func::Functional{:gga}, ρ::AbstractMatrix{T},
                      σ::AbstractMatrix{U}) where {T,U}
    @assert has_energy(func)  # Otherwise custom implementation of this function needed
    s_ρ, n_p = size(ρ)
    s_σ = size(σ, 1)
    TT = arithmetic_type(func, T, U)

    e   = similar(ρ, TT, n_p)
    Vρ  = similar(ρ, TT, s_ρ, n_p)
    Vσ  = similar(ρ, TT, s_σ, n_p)
    Vρρ = similar(ρ, TT, s_ρ, s_ρ, n_p)
    Vρσ = similar(ρ, TT, s_ρ, s_σ, n_p)
    Vσσ = similar(ρ, TT, s_σ, s_σ, n_p)

    # TODO Needed to make forward-diff work with !isbits floating-point types (e.g. BigFloat)
    Vρ  .= zero(TT)
    Vσ  .= zero(TT)
    Vρρ .= zero(TT)
    Vρσ .= zero(TT)
    Vσσ .= zero(TT)

    @views for i = 1:n_p
        kernel_terms!(e[i:i], Vρ[:, i], Vσ[:, i],
                      Vρρ[:, :, i], Vρσ[:, :, i], Vσσ[:, :, i],
                      func, ρ[:, i], σ[:, i])
    end
    (; e, Vρ, Vσ, Vρρ, Vρσ, Vσσ)
end
function kernel_terms!(e, Vρ, Vσ, Vρρ, Vρσ, Vσσ, func::Functional{:gga},
                       ρ::AbstractVector, σ::AbstractVector)
    res = ForwardDiff.hessian!(DiffResults.DiffResult(zero(eltype(e)), Vρ, Vρρ),
                               ρ -> energy(func, ρ, σ), ρ)
    res = ForwardDiff.hessian!(DiffResults.DiffResult(zero(eltype(e)), Vσ, Vσσ),
                               σ -> energy(func, ρ, σ), σ)

    dedρ = σ -> ForwardDiff.gradient(ρ -> energy(func, ρ, σ), ρ)
    ForwardDiff.jacobian!(Vρσ, dedρ, σ)

    e .= DiffResults.value(res)
    nothing
end

function energy(func::Functional{:gga}, ρ::AbstractVector{T},
                σ::AbstractVector{U}) where {T,U}
    length(ρ) == 1 || error("Multiple spins not yet implemented for fallback functionals")
    @assert length(ρ) == 1

    ρtotal = ρ[1]
    σtotal = σ[1]
    if ρtotal ≤ threshold_ρ(func, T)
        zero(arithmetic_type(func, T, U))
    else
        σstable = max(σtotal, threshold_σ(func, U))
        energy(func, ρtotal, σstable)
    end
end
