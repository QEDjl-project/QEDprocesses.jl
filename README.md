# QEDprocesses

[![Doc Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qedjl-project.github.io/QEDprocesses.jl/stable)
[![Doc Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qedjl-project.github.io/QEDprocesses.jl/dev)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Installation

To install the current stable version of `QEDprocesses.jl` you may use the standard julia package manager within the julia REPL

```julia
julia> using Pkg

julia> Pkg.add("QEDprocesses")
```

or you use the Pkg prompt by hitting `]` within the Julia REPL and then type

```julia
(@v1.9) pkg> add QEDprocesses
```

To install the locally downloaded package on Windows, change to the parent directory and type within the Pkg prompt

```julia
(@v1.9) pkg> add ./QEDprocesses.jl
```

## Building the documentation locally

To build the documentation of `QEDprocesses.jl` locally, first clone this
repository. Then, you instantiate the documentation subpackage by hitting 

```julia 
julia --project=docs -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
```
in the root directory of this repository. Afterwards, the dokumentation can be
built by running

```julia
julia --project=docs --color=yes docs/make.jl
```

To access the documentation site, just open the file `docs/_build/index.html` in
your favorite browser. 
