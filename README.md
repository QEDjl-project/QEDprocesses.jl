# QEDprocesses.jl

[![Doc Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qedjl-project.github.io/QEDprocesses.jl/stable)
[![Doc Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qedjl-project.github.io/QEDprocesses.jl/dev)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

**QEDprocesses.jl** is a Julia package for modeling scattering processes in quantum electrodynamics (QED).
It is part of the [`QuantumElectrodynamics.jl`](https://qedjl-project.github.io/QuantumElectrodynamics.jl/dev/) ecosystem,
which provides a modular framework for simulating QED interactions.

For details on how this package integrates with the rest of the `QuantumElectrodynamics.jl` ecosystem,
see the [QuantumElectrodynamics.jl documentation](https://qedjl-project.github.io/QuantumElectrodynamics.jl/dev/).

## Installation

To install the latest stable release of `QEDprocesses.jl`, use Julia's package manager:

```julia
using Pkg
Pkg.add("QEDprocesses")
```

Alternatively, in the Pkg prompt (accessed by pressing `]` in the Julia REPL):

```julia
(@v1.10) pkg> add QEDprocesses
```

To install a locally downloaded version of the package (e.g., on Windows), navigate to the parent directory and run:

```julia
(@v1.10) pkg> add ./QEDprocesses.jl
```

## Contributing

Contributions are welcome! If you find a bug, have a feature request, or want to contribute code, please open an issue or submit a pull request.

To maintain consistency across the `QuantumElectrodynamics.jl` ecosystem, contributors are encouraged to review the [QuantumElectrodynamics.jl development guide](https://qedjl-project.github.io/QuantumElectrodynamics.jl/stable/dev_guide/#Development-Guide).

## Credits and Contributors

This work was partly funded by the Center for Advanced Systems Understanding (CASUS) that
is financed by Germany’s Federal Ministry of Education and Research (BMBF) and by the Saxon
Ministry for Science, Culture and Tourism (SMWK) with tax funds on the basis of the budget
approved by the Saxon State Parliament.

### Core Contributors

- **Uwe Hernandez Acosta** (CASUS/HZDR, [u.hernandez@hzdr.de](mailto:u.hernandez@hzdr.de))
- **Anton Reinhard** (CASUS/HZDR)
- **Simeon Ehrig** (CASUS/HZDR)
- **Klaus Steiniger** (CASUS/HZDR)

### Former Contributors

- **Tom Jungnickel**

We sincerely thank all contributors who have supported this project.

### Acknowledgements

Special thanks to the following individuals for their support and contributions:

- **Michael Bussmann**
- **Tobias Dornheim**

## License

`QEDprocesses.jl` is licensed under the [MIT License](LICENSE).
