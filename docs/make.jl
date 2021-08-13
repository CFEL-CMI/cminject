using CMInject
using Documenter

DocMeta.setdocmeta!(CMInject, :DocTestSetup, :(using CMInject); recursive=true)

makedocs(;
    modules=[CMInject],
    authors="Simon Welker <simon.welker@cfel.de> and contributors",
    repo="https://github.com/CFEL-CMI/CMInject.jl/blob/{commit}{path}#{line}",
    sitename="CMInject.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://CFEL-CMI.github.io/CMInject.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/CFEL-CMI/CMInject.jl",
)
