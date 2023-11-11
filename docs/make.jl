using DftFunctionals
using Documenter

DocMeta.setdocmeta!(DftFunctionals, :DocTestSetup, :(using DftFunctionals); recursive=true)

makedocs(;
    modules=[DftFunctionals],
    authors="Michael F. Herbst <info@michael-herbst.com>",
    sitename="DftFunctionals",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliamolsim.github.io/DftFunctionals.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "interface.md",
        "generic.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/JuliaMolSim/DftFunctionals.jl",
    devbranch="master",
)
