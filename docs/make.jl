using DftFunctionals
using Documenter

DocMeta.setdocmeta!(DftFunctionals, :DocTestSetup, :(using DftFunctionals); recursive=true)

makedocs(;
    modules=[DftFunctionals],
    authors="Michael F. Herbst <info@michael-herbst.com>",
    repo="https://github.com/JuliaMolSim/DftFunctionals.jl/blob/{commit}{path}#{line}",
    sitename="DftFunctionals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliamolsim.github.io/DftFunctionals.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "interface.md",
        "generic.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaMolSim/DftFunctionals.jl",
    devbranch="master",
)
