struct LypCorrelation{CA} <: Functional{:gga,:c} where {CA<:ComponentArray{<:Number}}
    parameters::CA
    identifier::Symbol
end
function LypCorrelation(parameters::ComponentArray)
    LypCorrelation(parameters, :gga_c_lyp_custom)
end

identifier(pbe::LypCorrelation) = pbe.identifier
parameters(pbe::LypCorrelation) = pbe.parameters
function change_parameters(pbe::LypCorrelation, parameters::ComponentArray;
                           keep_identifier=false)
    if keep_identifier
        PbeCorrelation(parameters, pbe.identifier)
    else
        PbeCorrelation(parameters)
    end
end

function energy(lyp::LypCorrelation, ρ::T, σ::U) where {T<:Number,U<:Number}
    TT = arithmetic_type(lyp, T, U)

    # TODO This function is quite sensitive to the floating-point type ...
    #      so for now we don't bother doing this in TT, but rather convert before return
    #
    # TODO Comparing to libxc, this implementation has quite some numerical stability issues
    #      in the kernel
    #
    a = lyp.parameters.a
    b = lyp.parameters.b
    c = lyp.parameters.c
    d = lyp.parameters.d

    rr = cbrt(1 / ρ)  # (Wigner radius rₛ without prefactors)

    # Follows the partial integration trick in DOI 10.1016/0009-2614(89)87234-3, equation (2)
    δ(rr) = c * rr + d * rr / (1 + d * rr)
    ω(rr) =    exp(-c * rr) / (1 + d * rr)  # Factor ρ^(-11/3) absorbed into other terms
                                            # for stability reasons
    CF = 3/10 * (3*π^2)^(2/3)

    # ζ spin polarization (== 0 for non-spin-polarised)
    ζ = false

    # Again assuming no spin polarisation |∇ρ_{α}|² = |½ ∇ρ|² = ¼ |∇ρ|² = σ/4
    ∇ρα² = σ/4
    ∇ρβ² = σ/4
    ρα_ρ = (1 + ζ)/2  # = ρ_α / ρ
    ρβ_ρ = (1 - ζ)/2  # = ρ_β / ρ

    # Term in the square brackets in (2)
    term_square_bracket = (
         2.0^(11/3) * CF * ((ρα_ρ*ρ)^(8/3) + (ρβ_ρ*ρ)^(8/3))
       + (47/18 - 7/18 * δ(rr)) * σ
       - (5/2 - δ(rr)/18) * (∇ρα² + ∇ρβ²)
       - (δ(rr) - 11)/9 * ((1+ζ)/2 * ∇ρβ² + (1-ζ)/2 * ∇ρα²)
       #                    ρ_α/ρ            ρ_β/ρ
    )

    # Note, that (1 - ζ^2) = 4ρ_α ρ_β / ρ^2
    res = (
        -a * (1-ζ^2)*ρ / (1 + d*rr)
        -a*b * ω(rr) * rr^5 * (  term_square_bracket * (1 - ζ^2)/4
                               -  (2/3) * σ
                               + ((2/3) - ρα_ρ^2) * ∇ρβ²
                               + ((2/3) - ρβ_ρ^2) * ∇ρα² )
    )

    TT(res)
end

#
# Concrete functionals
#

"""
Standard LYP correlation.
Lee, Yang Parr 1988 (DOI: 10.1103/PhysRevB.37.785)
"""
function DftFunctional(::Val{:gga_c_lyp})
    LypCorrelation(ComponentArray(; a=0.04918, b=0.132, c=0.2533,  d=0.349), :gga_c_lyp)
end
