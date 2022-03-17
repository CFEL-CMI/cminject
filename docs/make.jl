using CMInject
using Documenter

DocMeta.setdocmeta!(CMInject, :DocTestSetup, :(using CMInject); recursive=true)

makedocs(;
    modules=[CMInject],
    authors="Simon Welker <simon.welker@cfel.de>, Timo Borner <timo.borner@cfel.de> and contributors",
    repo="https://github.com/CFEL-CMI/cminject/blob/{commit}{path}#{line}",
    sitename="CMInject.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://CFEL-CMI.github.io/cminject",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => "userguide.md",
        "Developer Guide" => "developerguide.md",
        "APIdocs" => "apidocs.md",
    ],
)

deploydocs(;
    repo="github.com/CFEL-CMI/cminject",
)
