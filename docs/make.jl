using QEDprocesses
using Documenter

DocMeta.setdocmeta!(QEDprocesses, :DocTestSetup, :(using QEDprocesses); recursive = true)

makedocs(;
    modules = [QEDprocesses],
    authors = "Uwe Hernandez Acosta <u.hernandez@hzdr.de>, Simeon Ehrig, Klaus Steiniger, Tom Jungnickel, Anton Reinhard",
    repo = Documenter.Remotes.GitHub("QEDjl-project", "QEDprocesses.jl"),
    sitename = "QEDprocesses.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://qedjl-project.gitlab.io/QEDprocesses.jl",
        edit_link = "dev",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)
deploydocs(repo = "github.com/QEDjl-project/QEDprocesses.jl.git", push_preview = false)
