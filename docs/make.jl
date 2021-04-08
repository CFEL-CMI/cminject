using CMInject
using Documenter

DocMeta.setdocmeta!(CMInject, :DocTestSetup, :(using CMInject); recursive=true)

makedocs(;
    modules=[CMInject],
    authors="Simon Welker <simon.welker@outlook.com> and contributors",
    repo="https://github.com/cobalamin/CMInject.jl/blob/{commit}{path}#{line}",
    sitename="CMInject.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cobalamin.github.io/CMInject.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cobalamin/CMInject.jl",
)
