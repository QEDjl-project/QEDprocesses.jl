using QEDprocesses
using Documenter

DocMeta.setdocmeta!(QEDprocesses, :DocTestSetup, :(using QEDprocesses); recursive=true)

makedocs(;
    modules=[QEDprocesses],
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de>, Simeon Ehrig, Klaus Steiniger, Tom Jungnickel, Anton Reinhard",
    repo="https://github.com/QEDjl-project/QEDprocesses.jl/blob/{commit}{path}#{line}",
    sitename="QEDprocesses.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
