using Dendrometry
using Documenter

DocMeta.setdocmeta!(Dendrometry, :DocTestSetup, :(using Dendrometry); recursive=true)

makedocs(;
    modules=[Dendrometry],
    authors="marcosdanieldasilva <marcosdasilva@5a.tec.br> and contributors",
    sitename="Dendrometry.jl",
    format=Documenter.HTML(;
        canonical="https://marcosdanieldasilva.github.io/Dendrometry.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/marcosdanieldasilva/Dendrometry.jl",
    devbranch="master",
)
