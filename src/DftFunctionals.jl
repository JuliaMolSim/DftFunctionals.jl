module DftFunctionals
using ForwardDiff

include("interface.jl")
include("DftFunctional.jl")
include("functionals/lda.jl")
include("functionals/gga_x_pbe.jl")
include("functionals/gga_c_pbe.jl")

export Functional, DftFunctional
export family, kind, identifier
export needs_σ, needs_τ, needs_Δρ, has_energy
export potential_terms, kernel_terms

export LdaExchange, PbeExchange, PbeCorrelation
end
