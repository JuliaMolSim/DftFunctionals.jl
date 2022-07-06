using Libxc_jll
const libxc = Libxc_jll.libxc

# Helper function to setup functionals in Libxc_jll
function xc_functional_alloc(identifier::Symbol; n_spin=1)
    pointer = ccall((:xc_func_alloc, libxc), Ptr{Cvoid}, ())
    number  = ccall((:xc_functional_get_number, libxc), Cint, (Cstring,), string(identifier))
    @assert number ≥ 0
    ret = ccall((:xc_func_init, libxc), Cint, (Ptr{Cvoid}, Cint, Cint), pointer, number, n_spin)
    @assert ret == 0
    pointer
end

function xc_functional_free(pointer)
    ccall((:xc_func_end, libxc), Cvoid, (Ptr{Cvoid},), pointer)
    ccall((:xc_func_free, libxc), Cvoid, (Ptr{Cvoid},), pointer)
end


function xc_lda(p, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4)
    ccall((:xc_lda, libxc), Cvoid, (Ptr{Cvoid}, Csize_t, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}),
          p, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4)
end

function xc_gga(p, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
                v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4)
    ccall((:xc_gga, libxc), Cvoid, (Ptr{Cvoid}, Csize_t, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                    Ptr{Cdouble}),
          p, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
          v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3, v4rho4, v4rho3sigma,
          v4rho2sigma2, v4rhosigma3, v4sigma4)
end
