const LdaExchange = DftFunctional{:lda_x,:lda,:x}

"""
LDA Slater exchange (DOI: 10.1017/S0305004100016108 and 10.1007/BF01340281)
"""
function energy(::LdaExchange, ρ::T) where {T <: Number}
    # Severe numerical issues if this is not done at least at Float64
    W = promote_type(Float64, T)
    T(-3/W(4) * cbrt(3/W(π) * ρ) * ρ)
end


# Wrapper for LDA correlation to correlation per particle (which is a primitive
# in GGA correlation functionals we want to keep numerical precision)
function energy(fun::Functional{:lda,:c}, ρ::T; kwargs...) where {T <: Number}
    energy_per_particle(fun, ρ; kwargs...) * ρ
end


"""
VWN5 LDA correlation according to Vosko, Wilk, and Nusair, (DOI 10.1139/p80-159).
"""
function energy_per_particle(::DftFunctional{:lda_c_vwn}, ρ::T) where {T <: Number}
    # From https://math.nist.gov/DFTdata/atomdata/node5.html
    A   = T( 0.0310907)
    x0  = T(-0.10498)
    b   = T( 3.72744)
    c   = T( 12.9352)
    rₛ  = cbrt(3/(T(4π)*ρ)) # τ in the above link
    x   = sqrt(rₛ)
    Xx  = x^2 + b*x + c
    Xx0 = x0^2 + b*x0 + c
    Q   = sqrt(4c-b^2)
    A * (log(x^2 / Xx) + 2b/Q*atan(Q/(2x+b)) - b*x0/Xx0*(log((x-x0)^2/Xx) + 2*(b+2x0)/Q*atan(Q/(2x+b))))
end



#
# LdaCorrelationPw
#
struct LdaCorrelationPw{Improved} <: Functional{:lda, :c}
end
LdaCorrelationPw(; improved=true) = LdaCorrelationPw{improved}()
identifier(pbe::LdaCorrelationPw) = :lda_c_pw
DftFunctional{:lda_c_pw}(;kwargs...) = LdaCorrelationPw(;kwargs...)


"""
Perdew, Wang correlation from 1992 (10.1103/PhysRevB.45.13244)
"""
function energy_per_particle(::LdaCorrelationPw{Improved}, ρ::T) where {T<:Number, Improved}
    α₁ = T.((0.21370,  0.20548,  0.11125))
    β₁ = T.((7.5957,  14.1189,  10.357  ))
    β₂ = T.((3.5876,   6.1977,   3.6231 ))
    β₃ = T.((1.6382,   3.3662,   0.88026))
    β₄ = T.((0.49294,  0.62517,  0.49671))
    # constant p = 1 is hard-coded in the expression for G below

    if !Improved  # Constants as given in original publication
        A   = T.((0.031091, 0.015545, 0.016887))
        f′′ = T(1.709921)  # f′′(0)
    else  # Modified constants, computed at improved accuracy
        A   = T.((0.0310907, 0.01554535, 0.0168869))
        f′′ = 8 / (9 * 2(cbrt(T(2)) - 1))  # f′′(0)
    end

    energy_per_particle_c_pw(ρ; A, α₁, β₁, β₂, β₃, β₄, f′′)
end
function energy_per_particle_c_pw(ρ::T; A, α₁, β₁, β₂, β₃, β₄, f′′) where {T}
    function G(sqrt_rₛ, A, α₁, β₁, β₂, β₃, β₄)  # (10) with p = 1 hard-coded
        denom = β₁ * sqrt_rₛ + β₂ * sqrt_rₛ^2 + β₃ * sqrt_rₛ^3 + β₄ * sqrt_rₛ^4
        -2A * (1 + α₁*sqrt_rₛ^2) * log(1 + 1 / (2A * denom) )
    end

    # equation (9)
    f(ζ) = ((1+ζ)^(4/T(3)) + (1-ζ)^(4/T(3)) - 2)/(2^(4/T(3)) - 2)  # == 0 for non-spin-polarised

    ε_0(rₛ) =  G(sqrt(rₛ), A[1], α₁[1], β₁[1], β₂[1], β₃[1], β₄[1])  # ε_c(rₛ, 0)
    ε_1(rₛ) =  G(sqrt(rₛ), A[2], α₁[2], β₁[2], β₂[2], β₃[2], β₄[2])  # ε_c(rₛ, 1)
    α(rₛ)   = -G(sqrt(rₛ), A[3], α₁[3], β₁[3], β₂[3], β₃[3], β₄[3])  # α_c(rₛ)

    # equation (8)
    ε(rₛ, ζ) = ε_0(rₛ) + α(rₛ) * f(ζ)/f′′ * (1 - ζ^4) + (ε_1(rₛ) - ε_0(rₛ)) * f(ζ) * ζ^4

    rₛ = cbrt(3 / (4T(π)  * ρ))  # equation (1)
    ε(rₛ, #= ζ = =# zero(T))
end
