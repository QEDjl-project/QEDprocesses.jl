using Pkg

# targeting the correct source code
# this asumes the make.jl script is located in QEDprocesses.jl/docs
project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")
Pkg.develop(; path = project_path)
# temporarily necessary because processes used to have a compat that is gone after the `develop` above
Pkg.update()

using Documenter
using QEDprocesses

DocMeta.setdocmeta!(QEDprocesses, :DocTestSetup, :(using QEDprocesses); recursive = true)

readme_path = joinpath(project_path, "README.md")
index_path = joinpath(project_path, "docs/src/index.md")
license_path = "https://github.com/QEDjl-project/QEDprocesses.jl/blob/dev/LICENSE"

# Copy README.md from the project base folder and use it as the start page
open(readme_path, "r") do readme_in
    readme_string = read(readme_in, String)

    # replace relative links in the README.md
    readme_string = replace(readme_string, "[MIT](LICENSE)" => "[MIT]($(license_path))")

    open(index_path, "w") do readme_out
        write(readme_out, readme_string)
    end
end

pages = [
    "Home" => "index.md",
    "API Reference" => "api.md",
]

try
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
        pages = pages,
    )
finally
    # doing some garbage collection
    @info "GarbageCollection: remove generated landing page"
    rm(index_path)
end

deploydocs(; repo = "github.com/QEDjl-project/QEDprocesses.jl.git", push_preview = false)
