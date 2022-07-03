"""PBE exchange."""
struct PbeExchange{Tk,Tu} <: Functional{:gga, :x} where {Tk,Tu}
    identifier::Symbol
    κ::Tk
    μ::Tu
    β::Tu  # Not used, only stored for convenience
end
function PbeExchange(identifier=:gga_x_pbe_custom; κ, β=nothing, μ=(β/3 * π^2))
    β = something(β, 3μ / π^2)
    PbeExchange(identifier, κ, μ, β)
end
identifier(pbe::PbeExchange) = pbe.identifier


function energy(pbe::PbeExchange, ρ::T, σ::U) where {T <: Number, U <: Number}
    TT = promote_type(T, U)
    κ = TT(pbe.κ)
    μ = TT(pbe.μ)

    pbe_x_f(s²) = 1 + κ - κ^2 / (κ + μ * s²)   # (14)
    # rₛ = cbrt(3 / (4π  * ρ))                 # page 2, left column, top
    # kF = cbrt(3π^2 * ρ)                      # page 2, left column, top
    # s  = sqrt(σ) / (2kF * ρ)                 # below (9)
    s² = σ / ( ρ^(4/3) * 2cbrt(3π^2) )^2

    energy(LdaExchange(), ρ) * pbe_x_f(s²)     # (10)
end


#
# Concrete functionals
#

# TODO Issues with the docstrings ...
# """
# Standard PBE exchange.
# Perdew, Burke, Ernzerhof 1996 (DOI: 10.1103/PhysRevLett.77.3865)
# """
DftFunctional{:gga_x_pbe}() = PbeExchange(:gga_x_pbe; κ=0.8040, β=0.06672455060314922)

# """
# Revised PBE exchange.
# Zhang, Yang 1998 (DOI 10.1103/physrevlett.80.890)
# """
DftFunctional{:gga_x_pbe_r}() = PbeExchange(:gga_x_pbe_r; κ=1.245, β=0.06672455060314922)

# """
# XPBE exchange.
# Xu, Goddard 2004 (DOI 10.1063/1.1771632)
# """
DftFunctional{:gga_x_xpbe}() = PbeExchange(:gga_x_xpbe; κ=0.91954, μ=0.23214)  # Table 1

# """
# PBESol exchange.
# Perdew, Ruzsinszky, Csonka and others 2008 (DOI 10.1103/physrevlett.100.136406)
# """
DftFunctional{:gga_x_pbe_sol}() = PbeExchange(:gga_x_pbe_sol; κ=0.8040, μ=10/81)
# μ given below equation (2)

# """
# APBE exchange.
# Constantin, Fabiano, Laricchia 2011 (DOI 10.1103/physrevlett.106.186406)
# """
DftFunctional{:gga_x_apbe}() = PbeExchange(:gga_x_apbe; κ=0.8040, μ=0.260)
# p. 1, right column, bottom

# """
# PBEmol exchange.
# del Campo, Gazqez, Trickey and others 2012 (DOI 10.1063/1.3691197)
# """
DftFunctional{:gga_x_pbe_mol}() = PbeExchange(:gga_x_pbe_mol; κ=0.8040, μ=0.27583)
# p. 4, left column, bottom

# """
# PBEfe exchange.
# Sarmiento-Perez, Silvana, Marques 2015 (DOI 10.1021/acs.jctc.5b00529)
# """
DftFunctional{:gga_x_pbefe}() = PbeExchange(:gga_x_pbefe; κ=0.437, μ=0.346)  # Table 1
