# Changelog

## Version 0.2.0

[diff since 0.1.0](https://github.com/QEDjl-project/QEDprocesses.jl/compare/v0.1.0...release-0.2.0)

This release is part of the restructuring processes of QED.jl (see https://github.com/QEDjl-project/QED.jl/issues/35 for details). 
It is a breaking release, indicated by the bumped minor version because we don't have a major version for this yet.

### Breaking changes

This release moves the general purpose interfaces to [`QEDbase.jl`](https://github.com/QEDjl-project/QEDbase.jl) and the core functionality to the new package [`QEDcore.jl`](https://github.com/QEDjl-project/QEDcore.jl).
The purpose of this package is to provide implementations for specific models and processes, currently the 1-photon-Compton in perturbative QED.

### New features

- [#28](https://github.com/QEDjl-project/QEDprocesses.jl/pull/28): Add particle propagators.
- [#39](https://github.com/QEDjl-project/QEDprocesses.jl/pull/39): Coordinate-based implementation of diff cs and probabilities. (Superseded by later `PhaseSpacePoint` implementations, which can be constructed using coordinates)
- [#40](https://github.com/QEDjl-project/QEDprocesses.jl/pull/40): Add perturbative one-photon Compton implementation.
- [#51](https://github.com/QEDjl-project/QEDprocesses.jl/pull/51): Add initial `PhaseSpacePoint` definition.
- [#52](https://github.com/QEDjl-project/QEDprocesses.jl/pull/52): Add overloads in the process interface using `PhaseSpacePoint`s.
- [#54](https://github.com/QEDjl-project/QEDprocesses.jl/pull/54): Add phase space generation interface.
- [#59](https://github.com/QEDjl-project/QEDprocesses.jl/pull/59): Update process interface to use `PhaseSpacePoint`s instead of raw momenta. (Supersedes [#52](https://github.com/QEDjl-project/QEDprocesses.jl/pull/52))
- [#60](https://github.com/QEDjl-project/QEDprocesses.jl/pull/60): Add pretty-printing for some of the types. (`ParticleStateful`, `PhaseSpacePoint`, models, processes, coordinate systems, frames of reference, and phase space definitions)
- [#63](https://github.com/QEDjl-project/QEDprocesses.jl/pull/63): Improvement for the `PhaseSpacePoint` type, making it type stable and allowing in-place construction.
- [#68](https://github.com/QEDjl-project/QEDprocesses.jl/pull/68): Add an `InPhaseSpacePoint` definition allowing a phase space point to not contain particle states for the outgoing particles. This is useful in some interfaces that don't need the outgoing part of a phase space point.

### Maintenance

Besides the new features, this release contains some maintenance and minor changes and fixes. Many of these were temporarily necessary to maintain a working state during the restructuring.

- [#36](https://github.com/QEDjl-project/QEDprocesses.jl/pull/36): Remove custom registry from unit tests.
- [#37](https://github.com/QEDjl-project/QEDprocesses.jl/pull/37): Remove custom registry build step.
- [#38](https://github.com/QEDjl-project/QEDprocesses.jl/pull/38): Refactor differential cross section functionality.
- [#53](https://github.com/QEDjl-project/QEDprocesses.jl/pull/53): CompatHelper: StaticArrays compat version 1.
- [#61](https://github.com/QEDjl-project/QEDprocesses.jl/pull/61): CompatHelper: QuadGK compat version 2.
- [#62](https://github.com/QEDjl-project/QEDprocesses.jl/pull/62): Move test dependencies to the main `Project.toml`.
- [#71](https://github.com/QEDjl-project/QEDprocesses.jl/pull/71): Use Julia 1.10 in CI.
- [#74](https://github.com/QEDjl-project/QEDprocesses.jl/pull/74): (Re-)move interfaces that move to QEDbase.
- [#76](https://github.com/QEDjl-project/QEDprocesses.jl/pull/76): Refactor QED.jl package dependencies and namespaces.
- [#78](https://github.com/QEDjl-project/QEDprocesses.jl/pull/78): Use function implementations from QEDcore.
- [#80](https://github.com/QEDjl-project/QEDprocesses.jl/pull/80): Fix ambiguous function calls.
- [#81](https://github.com/QEDjl-project/QEDprocesses.jl/pull/81): Update dev version.
- [#82](https://github.com/QEDjl-project/QEDprocesses.jl/pull/82): Remove QEDbase compat entry.
- [#87](https://github.com/QEDjl-project/QEDprocesses.jl/pull/87): (Re-)move abstract cross section and probability implementations (move to QEDcore).
- [#89](https://github.com/QEDjl-project/QEDprocesses.jl/pull/89): Removes temporary dependencies to other dev versions that were necessary during restructuring.

## Version 0.1.0

This is the initial verison of QEDprocesses.

[Diff since inital commit](https://github.com/QEDjl-project/QEDprocesses.jl/compare/302274695d82225f4a810c252d6919839bc59fd7...release-v0.1.0) 
[Full list of PRs](https://github.com/QEDjl-project/QEDprocesses.jl/milestone/2?closed=1)


### Highlights
- interface for scattering processes described by physical models (general properties, differential/total cross sections) https://github.com/QEDjl-project/QEDprocesses.jl/pull/11
- interface for abstract setups (general computation setup and process dedicated
setups) https://github.com/QEDjl-project/QEDprocesses.jl/pull/14
- setup ci for unit and integration testing https://github.com/QEDjl-project/QEDprocesses.jl/pull/5 https://github.com/QEDjl-project/QEDprocesses.jl/pull/7 https://github.com/QEDjl-project/QEDprocesses.jl/pull/21
- setup ci/cd for docs https://github.com/QEDjl-project/QEDprocesses.jl/pull/19
