using DftFunctionals
using Documenter

DocMeta.setdocmeta!(DftFunctionals, :DocTestSetup, :(using DftFunctionals); recursive=true)

makedocs(;
    modules=[DftFunctionals],
    authors="Michael F. Herbst <info@michael-herbst.com> and contributors",
    repo="https://github.com/mfherbst/DftFunctionals.jl/blob/{commit}{path}#{line}",
    sitename="DftFunctionals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfherbst.github.io/DftFunctionals.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mfherbst/DftFunctionals.jl",
    devbranch="master",
)
