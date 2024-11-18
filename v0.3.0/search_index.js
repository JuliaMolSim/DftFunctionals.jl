var documenterSearchIndex = {"docs":
[{"location":"interface/","page":"Interface","title":"Interface","text":"CurrentModule = DftFunctionals","category":"page"},{"location":"interface/#Interface","page":"Interface","title":"Interface","text":"","category":"section"},{"location":"interface/","page":"Interface","title":"Interface","text":"API documentation for DftFunctionals.","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"See in particular the key functions potential_terms and kernel_terms.","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"note: Missing documentation\nMore details on the public interface are needed.","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"Modules = [DftFunctionals]\nPrivate = false","category":"page"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Symbol, Vararg{Any}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"DftFunctional(identifier::Symbol, args...; kwargs...)\n\nA generic DFT functional implementation. Valid identifiers are the ones supported by Libxc as symbols, e.g. :lda_x, hyb_gga_xc_b3lyp. Only a subset of functionals is currently available. Additional arguments and kwargs are passed to the functional (e.g. to modify functional parameters).\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_c_apbe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"APBE correlation. Constantin, Fabiano, Laricchia 2011 (DOI 10.1103/physrevlett.106.186406)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_c_pbe_mol}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"PBEmol correlation. del Campo, Gazqez, Trickey and others 2012 (DOI 10.1063/1.3691197)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_c_pbe_sol}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"PBESol correlation. Perdew, Ruzsinszky, Csonka and others 2008 (DOI 10.1103/physrevlett.100.136406)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_c_pbefe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"PBEfe correlation. Sarmiento-Perez, Silvana, Marques 2015 (DOI 10.1021/acs.jctc.5b00529)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_c_pbe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"Standard PBE correlation. Perdew, Burke, Ernzerhof 1996 (DOI: 10.1103/PhysRevLett.77.3865)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_c_xpbe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"XPBE correlation. Xu, Goddard 2004 (DOI 10.1063/1.1771632)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_apbe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"APBE exchange. Constantin, Fabiano, Laricchia 2011 (DOI 10.1103/physrevlett.106.186406)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_pbe_mol}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"PBEmol exchange. del Campo, Gazqez, Trickey and others 2012 (DOI 10.1063/1.3691197)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_pbe_r}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"Revised PBE exchange. Zhang, Yang 1998 (DOI 10.1103/physrevlett.80.890)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_pbe_sol}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"PBESol exchange. Perdew, Ruzsinszky, Csonka and others 2008 (DOI 10.1103/physrevlett.100.136406)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_pbefe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"PBEfe exchange. Sarmiento-Perez, Silvana, Marques 2015 (DOI 10.1021/acs.jctc.5b00529)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_pbe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"Standard PBE exchange. Perdew, Burke, Ernzerhof 1996 (DOI: 10.1103/PhysRevLett.77.3865)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:gga_x_xpbe}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"XPBE exchange. Xu, Goddard 2004 (DOI 10.1063/1.1771632)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:lda_c_pw}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"Perdew, Wang correlation from 1992 (10.1103/PhysRevB.45.13244)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:lda_c_vwn}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"VWN5 LDA correlation according to Vosko, Wilk, and Nusair, (DOI 10.1139/p80-159).\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.DftFunctional-Tuple{Val{:lda_x}}","page":"Interface","title":"DftFunctionals.DftFunctional","text":"LDA Slater exchange (DOI: 10.1017/S0305004100016108 and 10.1007/BF01340281)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.change_parameters","page":"Interface","title":"DftFunctionals.change_parameters","text":"Return a new version of the passed functional with its parameters adjusted. This may not be a copy in case no changes are done to its internal parameters. Generally the identifier of the functional will be changed to reflect the change in parameter values unless keep_identifier is true. To get the tuple of adjustable parameters and their current values check out parameters. It is not checked that the correct parameters are passed.\n\nchange_parameters(f::Functional, params_new; keep_identifier=false)::Functional\n\n\n\n\n\n","category":"function"},{"location":"interface/#DftFunctionals.family-Union{Tuple{Functional{F}}, Tuple{F}} where F","page":"Interface","title":"DftFunctionals.family","text":"Return the family of a functional. Results are :lda, :gga, :mgga and :mggal (Meta-GGA requiring Laplacian of ρ)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.has_energy-Tuple{Functional}","page":"Interface","title":"DftFunctionals.has_energy","text":"Does this functional support energy evaluations? Some don't, in which case energy terms will not be returned by potential_terms and kernel_terms, i.e. e will be false (a strong zero).\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.identifier","page":"Interface","title":"DftFunctionals.identifier","text":"Return the identifier corresponding to a functional\n\n\n\n\n\n","category":"function"},{"location":"interface/#DftFunctionals.kernel_terms","page":"Interface","title":"DftFunctionals.kernel_terms","text":"kernel_terms(f::Functional, ρ, [σ, τ, Δρ])\n\nEvaluate energy, potential and kernel terms at a real-space grid of densities, density derivatives etc. Not required derivatives for the functional type will be ignored. Returns a named tuple with the same keys as potential_terms and additionally second-derivative cross terms such as Vρσ (frac^2eρσ).\n\n\n\n\n\n","category":"function"},{"location":"interface/#DftFunctionals.kind-Union{Tuple{Functional{F, K}}, Tuple{K}, Tuple{F}} where {F, K}","page":"Interface","title":"DftFunctionals.kind","text":"Return the functional kind: :x (exchange), :c (correlation), :k (kinetic) or :xc (exchange and correlation combined)\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.needs_Δρ-Union{Tuple{Functional{F}}, Tuple{F}} where F","page":"Interface","title":"DftFunctionals.needs_Δρ","text":"True if the functional needs Δ ρ.\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.needs_σ-Union{Tuple{Functional{F}}, Tuple{F}} where F","page":"Interface","title":"DftFunctionals.needs_σ","text":"True if the functional needs σ = ρ  ρ.\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.needs_τ-Union{Tuple{Functional{F}}, Tuple{F}} where F","page":"Interface","title":"DftFunctionals.needs_τ","text":"True if the functional needs τ (kinetic energy density).\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.parameters-Tuple{Functional}","page":"Interface","title":"DftFunctionals.parameters","text":"Return adjustable parameters of the functional and their values.\n\n\n\n\n\n","category":"method"},{"location":"interface/#DftFunctionals.potential_terms","page":"Interface","title":"DftFunctionals.potential_terms","text":"potential_terms(f::Functional, ρ, [σ, τ, Δρ])\n\nEvaluate energy and potential terms at a real-space grid of densities, density derivatives etc. Not required derivatives for the functional type will be ignored. Returns a named tuple with keys e (Energy per unit volume), Vρ (fraceρ), Vσ (fraceσ), Vτ (fraceτ), Vl (frace(Δρ)).\n\n\n\n\n\n","category":"function"},{"location":"interface/#DftFunctionals.spinindex_σ-Tuple{Int64, Int64}","page":"Interface","title":"DftFunctionals.spinindex_σ","text":"In the functional interface we explicitly use the symmetry σ_αβ = σ_βα for the contracted density gradient σ_ij = ρ_i  ρ_j where ij  α β. Instead of treating σ to have the spin components σ_αα, σ_αβ, σ_βα and σ_ββ, we consider it to have three pseudo-spin components σ_αα, σ_x and σ_ββ, where σ_x = (σ_αβ + σ_βα)2 = σ_αβ. Input arrays like σ or output arrays like Vσ or Vρσ will therefore feature a spin axis of length 3 refering to the σ_αα, σ_αβ and σ_ββ components / derivatives respectively. E.g. σ is of shape (3, n_p) (where n_p is the number of points where evaluation takes place.\n\nThis function maps the \"unfolded\" spin tuple (i, j) in σ_ij to the corresponding place in the length-3 spin axis of σ.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DftFunctionals","category":"page"},{"location":"#DftFunctionals","page":"Home","title":"DftFunctionals","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DftFunctionals provides an Interface for evaluating exchange-correlation functionals for density-functional theory calculations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a small number of functionals a Generic functional implementation in pure Julia is also provided, where the focus for the moment is on readable code rather than speed.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package is at an early stage, so please report bugs or numerical issues if you find them.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A more stable (but less generic) implementation of DFT functionals is available via the Libxc.jl package. Currently Libxc is not yet compatible with DftFunctionals out of the box (but there exists an experimental wrapper in DFTK.jl). This will change in the future.","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Input data to DftFunctionals is the density ρ (and its derivatives such as σ = nabla rho cdot nabla rho etc.) on a real-space grid. This data needs to be passed as a (n_spin, n_p) matrix, where n_spin is the number of spin components and n_p is the number of data points. Notice that three-dimensional density therefore need to be flattened before calling anything from DftFunctionals.jl. For more details on the shapes of the input and output arrays, see the documentation of potential_terms.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using DftFunctionals\n\n# Setup dummy input data\nρ = reshape([0.1, 0.2, 0.3, 0.4, 0.5], 1, :)\nσ = reshape([0.2, 0.3, 0.4, 0.5, 0.6], 1, :)\n\n# Setup the functional to use\nxc = DftFunctional(:gga_x_pbe)\n\n# Get the potential terms\nterms = potential_terms(xc, ρ, σ)\n\n# Energy per unit volume on each grid point\n@show terms.e\n\n# Potential Vσ = de / dρ on each grid point\n@show terms.Vρ\n\nnothing","category":"page"},{"location":"generic/#Generic-functional-implementation","page":"Generic functional implementation","title":"Generic functional implementation","text":"","category":"section"},{"location":"generic/","page":"Generic functional implementation","title":"Generic functional implementation","text":"Generic implementations are available for the following functional identifier taken from Libxc. You can use them by constructing DftFunctional(<identifier>).","category":"page"},{"location":"generic/","page":"Generic functional implementation","title":"Generic functional implementation","text":":lda_x\n:lda_c_vwn\n:lda_c_pw\n:gga_x_pbe\n:gga_x_pbe_r\n:gga_x_xpbe\n:gga_x_pbe_sol\n:gga_x_apbe\n:gga_x_pbe_mol\n:gga_x_pbefe\n:gga_c_pbe\n:gga_c_xpbe\n:gga_c_pbe_sol\n:gga_c_apbe\n:gga_c_pbe_mol\n:gga_c_pbefe","category":"page"},{"location":"generic/","page":"Generic functional implementation","title":"Generic functional implementation","text":"note: Improve this\nLink to DOI or something or best generate this from the source directly.","category":"page"},{"location":"generic/#Implementing-new-functionals","page":"Generic functional implementation","title":"Implementing new functionals","text":"","category":"section"},{"location":"generic/","page":"Generic functional implementation","title":"Generic functional implementation","text":"note: Missing documentation\nMore details on how to implement new functionals.","category":"page"}]
}