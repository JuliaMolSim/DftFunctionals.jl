struct PbeExchange{Tκ, Tμ} <: Functional{:gga,:x} where {Tκ<:Number,Tμ<:Number}
    κ::Tκ
    μ::Tμ
end
function PbeExchange(; κ, μ)
    PbeExchange(κ, μ)
end
function parameters_type(pbe::PbeExchange)
    promote_type(typeof(pbe.κ), typeof(pbe.μ))
end

function energy(pbe::PbeExchange, ρ::T, σ::U) where {T<:Number,U<:Number}
    TT = arithmetic_type(pbe, T, U)

    # TODO This function is quite sensitive to the floating-point type ...
    #      so for now we don't bother doing this in TT, but rather convert before return
    κ = pbe.κ
    μ = pbe.μ

    pbe_x_f(s²) = 1 + κ - κ^2 / (κ + μ * s²)   # (14)
    # rₛ = cbrt(3 / (4π  * ρ))                 # page 2, left column, top
    # kF = cbrt(3π^2 * ρ)                      # page 2, left column, top
    # s  = sqrt(σ) / (2kF * ρ)                 # below (9)
    s² = σ / (ρ^(4 / 3) * 2cbrt(3π^2))^2

    res = energy(LdaExchange(), ρ) * pbe_x_f(s²)     # (10)
    TT(res)
end

# Conversion between μ and β (some authors use one, some the other)
pbe_μ_from_β(β) = β / 3 * π^2
pbe_β_from_μ(μ) = 3μ / π^2

#
# Concrete functionals
#

const KNOWN_X_PBE = [
    # Standard PBE exchange.
    # Perdew, Burke, Ernzerhof 1996 (DOI: 10.1103/PhysRevLett.77.3865)
    :gga_x_pbe     => (; κ=0.8040, μ=pbe_μ_from_β(0.06672455060314922)),
    # Revised PBE exchange.
    # Zhang, Yang 1998 (DOI 10.1103/physrevlett.80.890)
    :gga_x_pbe_r   => (; κ=1.245, μ=pbe_μ_from_β(0.06672455060314922)),
    # XPBE exchange.
    # Xu, Goddard 2004 (DOI 10.1063/1.1771632)
    :gga_x_xpbe    => (; κ=0.91954, μ=0.23214), # Table 1
    # PBESol exchange.
    # Perdew, Ruzsinszky, Csonka and others 2008 (DOI 10.1103/physrevlett.100.136406)
    # μ given below equation (2)
    :gga_x_pbe_sol => (; κ=0.8040, μ=10 / 81),
    # APBE exchange.
    # Constantin, Fabiano, Laricchia 2011 (DOI 10.1103/physrevlett.106.186406)
    # p. 1, right column, bottom
    :gga_x_apbe    => (; κ=0.8040, μ=0.260),
    # PBEmol exchange.
    # del Campo, Gazqez, Trickey and others 2012 (DOI 10.1063/1.3691197)
    # p. 4, left column, bottom
    :gga_x_pbe_mol => (; κ=0.8040, μ=0.27583),
    # PBEfe exchange.
    # Sarmiento-Perez, Silvana, Marques 2015 (DOI 10.1021/acs.jctc.5b00529)
    :gga_x_pbefe   => (; κ=0.437, μ=0.346), # Table 1
]

for (id, param) in KNOWN_X_PBE
    @eval function DftFunctional(::Val{$(QuoteNode(id))})
        PbeExchange(κ=$(param.κ), μ=$(param.μ))
    end
end

function identifier(pbe::PbeExchange)
    for (id, param) in KNOWN_X_PBE
        if pbe.κ ≈ param.κ && pbe.μ ≈ param.μ
            return id
        end
    end
    nothing
end
