module DftFunctionals
using ForwardDiff
using ComponentArrays

include("interface.jl")
include("util.jl")
export Functional
export family, kind, identifier
export parameters, change_parameters
export needs_σ, needs_τ, needs_Δρ, has_energy
export potential_terms, kernel_terms
export spinindex_σ

# Generic implementation
include("DftFunctional.jl")
include("functionals/lda.jl")
include("functionals/gga_x_pbe.jl")
include("functionals/gga_c_pbe.jl")
export DftFunctional
export LdaExchange, LdaCorrelationVwn, LdaCorrelationPw, PbeExchange, PbeCorrelation

end
