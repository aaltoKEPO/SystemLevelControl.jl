using SystemLevelControl
using Documenter

DocMeta.setdocmeta!(SystemLevelControl, :DocTestSetup, :(using SystemLevelControl); recursive=true)

makedocs(;
    modules=[SystemLevelControl],
    authors="Otacilio 'Minho' Neto <otacilio.neto.aalto.fi>",
    repo="https://github.com/tiominho/SystemLevelControl.jl/blob/{commit}{path}#{line}",
    sitename="SystemLevelControl.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tiominho.github.io/SystemLevelControl.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tiominho/SystemLevelControl.jl",
    devbranch="main",
)
