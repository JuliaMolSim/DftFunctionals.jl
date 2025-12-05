struct PbeCorrelation{Tlda,Tβ,Tγ} <:
       Functional{:gga,:c} where {Tlda,Tβ<:Number,Tγ<:Number}
    lda::Tlda
    β::Tβ
    γ::Tγ
end
function PbeCorrelation(; lda=DftFunctional(:lda_c_pw), β, γ)
    PbeCorrelation(lda, β, γ)
end

function parameters_type(pbe::PbeCorrelation)
    promote_type(parameters_type(pbe.lda), typeof(pbe.β), typeof(pbe.γ))
end

function energy(pbe::PbeCorrelation, ρ::T, σ::U) where {T<:Number,U<:Number}
    TT = arithmetic_type(pbe, T, U)

    # TODO This function is quite sensitive to the floating-point type ...
    #      so for now we don't bother doing this in TT, but rather convert before return
    β = pbe.β
    γ = pbe.γ

    # Spin-scaling factor with ζ spin polarization.
    # Yue Wang and John P. Perdew. Phys. Rev. B 43, 8911 (1991).
    # DOI 10.1103/PhysRevB.43.8911
    ϕ(ζ) = ((1 + ζ)^(2 / 3) + (1 - ζ)^(2 / 3)) / 2  # == 1 for non-spin-polarised

    # ε = UEG correlation energy per particle
    A(ε, ϕ³) = β / γ / expm1(-ε / (γ * ϕ³))  # (8)
    function H(ε, t², ϕ³)  # (7)
        At² = A(ε, ϕ³) * t²
        γ * ϕ³ * log(1 + β / γ * t² * (1 + At²) / (1 + At² + (At²)^2))
    end

    phi = 1.0 #= ϕ(ζ) =#
    ε_lda = energy_per_particle(pbe.lda, ρ)
    t² = (1 / 12 * 3^(5 / 6) * π^(1 / 6))^2 * σ / (phi^2 * ρ^(7 / 3))  # page 2, left column, top
    res = (ε_lda + H(ε_lda, t², phi^3)) * ρ

    TT(res)
end

#
# Concrete functionals
#

const KNOWN_C_PBE = [
    # Standard PBE correlation.
    # Perdew, Burke, Ernzerhof 1996 (DOI: 10.1103/PhysRevLett.77.3865)
    :gga_c_pbe => (; β=0.06672455060314922, γ=(1 - log(2)) / π^2),
    # XPBE correlation.
    # Xu, Goddard 2004 (DOI 10.1063/1.1771632)
    :gga_c_xpbe => let
        β = 0.089809  # Fitted constants, Table I
        α = 0.197363  # Fitted constants, Table I
        γ = β^2 / 2α
        (; β, γ)
    end,
    # PBESol correlation.
    # Perdew, Ruzsinszky, Csonka and others 2008 (DOI 10.1103/physrevlett.100.136406)
    # Page 3, left column below figure 1
    :gga_c_pbe_sol => (; β=0.046, γ=(1 - log(2)) / π^2),
    # APBE correlation.
    # Constantin, Fabiano, Laricchia 2011 (DOI 10.1103/physrevlett.106.186406)
    :gga_c_apbe => let
        μ = 0.260   # p. 1, right column, bottom
        β = 3μ / π^2
        γ = (1 - log(2)) / π^2  # like in PBE
        (; β, γ)
    end,
    # PBEmol correlation.
    # del Campo, Gazqez, Trickey and others 2012 (DOI 10.1063/1.3691197)
    # β made to cancel self-interaction error in hydrogen
    # p. 4, right column, first paragraph
    :gga_c_pbe_mol => (; β=0.08384, γ=(1 - log(2)) / π^2),
    # PBEfe correlation.
    # Sarmiento-Perez, Silvana, Marques 2015 (DOI 10.1021/acs.jctc.5b00529)
    # Fitted constants, Table I
    :gga_c_pbefe => (; β=0.043, γ=0.031090690869654895034),
]

for (id, param) in KNOWN_C_PBE
    @eval function DftFunctional(::Val{$(QuoteNode(id))})
        PbeCorrelation(β=$(param.β), γ=$(param.γ))
    end
end

function identifier(pbe::PbeCorrelation)
    for (id, param) in KNOWN_C_PBE
        if pbe.β ≈ param.β && pbe.γ ≈ param.γ
            return id
        end
    end
    nothing
end
