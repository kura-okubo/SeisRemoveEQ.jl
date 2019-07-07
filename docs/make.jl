using Documenter, SeisRemoveEQ

makedocs(;
    modules=[SeisRemoveEQ],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/kura-okubo/SeisRemoveEQ.jl/blob/{commit}{path}#L{line}",
    sitename="SeisRemoveEQ.jl",
    authors="kurama",
    assets=String[],
)

deploydocs(;
    repo="github.com/kura-okubo/SeisRemoveEQ.jl",
)
