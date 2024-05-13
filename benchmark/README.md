# Benchmarks for QEDprocesses.jl

This folder contains some micro benchmarks for the core functionality of
`QEDprocesses.jl`.

## Run benchmarks

To run the benchmarks locally on your machine, go to the subfolder `benchmark` and execute
the script `run_local.jl`:

```console
cd benchmark
julia run_local.jl
```

This builds the benchmark suite, tunes the benchmarks, runs them and saves the results
locally in a file `bench.json`.

## Plot results

For plotting the results, one can just run

```console
julia plot_bench.jl
```

which creates all benchmark plots and stores them in the folder `plots`.
