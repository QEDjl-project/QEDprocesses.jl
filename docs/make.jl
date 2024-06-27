using Pkg

# targeting the correct source code
# this asumes the make.jl script is located in QEDprocesses.jl/docs
project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")
Pkg.develop(; path=project_path)
# temporarily necessary because processes used to have a compat that is gone after the `develop` above
Pkg.update()

using Documenter
using QEDprocesses

DocMeta.setdocmeta!(QEDprocesses, :DocTestSetup, :(using QEDprocesses); recursive=true)

makedocs(;
    modules=[QEDprocesses],
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de>, Simeon Ehrig, Klaus Steiniger, Tom Jungnickel, Anton Reinhard",
    repo=Documenter.Remotes.GitHub("QEDjl-project", "QEDprocesses.jl"),
    sitename="QEDprocesses.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://qedjl-project.gitlab.io/QEDprocesses.jl",
        edit_link="dev",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)
deploydocs(; repo="github.com/QEDjl-project/QEDprocesses.jl.git", push_preview=false)
