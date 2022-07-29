struct PbeExchange{CA} <: Functional{:gga, :x} where {CA <: ComponentArray}
    identifier::Symbol
    parameters::CA
end
PbeExchange(parameters::ComponentArray) = PbeExchange(:gga_x_pbe_custom, parameters)

identifier(pbe::PbeExchange) = pbe.identifier
parameters(pbe::PbeExchange) = pbe.parameters
function change_parameters(pbe::PbeExchange, parameters::ComponentArray; keep_identifier=false)
    if keep_identifier
        PbeExchange(pbe.identifier, parameters)
    else
        PbeExchange(parameters)
    end
end


function energy(pbe::PbeExchange, ρ::T, σ::U) where {T <: Number, U <: Number}
    TT = promote_type(T, U, parameter_type(pbe))
    κ = TT(pbe.parameters.κ)
    μ = TT(pbe.parameters.μ)

    pbe_x_f(s²) = 1 + κ - κ^2 / (κ + μ * s²)   # (14)
    # rₛ = cbrt(3 / (4π  * ρ))                 # page 2, left column, top
    # kF = cbrt(3π^2 * ρ)                      # page 2, left column, top
    # s  = sqrt(σ) / (2kF * ρ)                 # below (9)
    s² = σ / ( ρ^(4/3) * 2cbrt(3π^2) )^2

    energy(LdaExchange(), ρ) * pbe_x_f(s²)     # (10)
end

# Conversion between μ and β (some authors use one, some the other)
pbe_μ_from_β(β) = β/3 * π^2
pbe_β_from_μ(μ) = 3μ / π^2


#
# Concrete functionals
#

"""
Standard PBE exchange.
Perdew, Burke, Ernzerhof 1996 (DOI: 10.1103/PhysRevLett.77.3865)
"""
function DftFunctional(::Val{:gga_x_pbe})
    PbeExchange(:gga_x_pbe, ComponentArray(κ=0.8040, μ=pbe_μ_from_β(0.06672455060314922)))
end

"""
Revised PBE exchange.
Zhang, Yang 1998 (DOI 10.1103/physrevlett.80.890)
"""
function DftFunctional(::Val{:gga_x_pbe_r})
    PbeExchange(:gga_x_pbe_r, ComponentArray(κ=1.245, μ=pbe_μ_from_β(0.06672455060314922)))
end

"""
XPBE exchange.
Xu, Goddard 2004 (DOI 10.1063/1.1771632)
"""
function DftFunctional(::Val{:gga_x_xpbe})
    PbeExchange(:gga_x_xpbe, ComponentArray(κ=0.91954, μ=0.23214))  # Table 1
end

"""
PBESol exchange.
Perdew, Ruzsinszky, Csonka and others 2008 (DOI 10.1103/physrevlett.100.136406)
"""
function DftFunctional(::Val{:gga_x_pbe_sol})
    # μ given below equation (2)
    PbeExchange(:gga_x_pbe_sol, ComponentArray(κ=0.8040, μ=10/81))
end

"""
APBE exchange.
Constantin, Fabiano, Laricchia 2011 (DOI 10.1103/physrevlett.106.186406)
"""
function DftFunctional(::Val{:gga_x_apbe})
    # p. 1, right column, bottom
    PbeExchange(:gga_x_apbe, ComponentArray(κ=0.8040, μ=0.260))
end

"""
PBEmol exchange.
del Campo, Gazqez, Trickey and others 2012 (DOI 10.1063/1.3691197)
"""
function DftFunctional(::Val{:gga_x_pbe_mol})
    # p. 4, left column, bottom
    PbeExchange(:gga_x_pbe_mol, ComponentArray(κ=0.8040, μ=0.27583))
end

"""
PBEfe exchange.
Sarmiento-Perez, Silvana, Marques 2015 (DOI 10.1021/acs.jctc.5b00529)
"""
function DftFunctional(::Val{:gga_x_pbefe})
    PbeExchange(:gga_x_pbefe, ComponentArray(κ=0.437, μ=0.346))  # Table 1
end
