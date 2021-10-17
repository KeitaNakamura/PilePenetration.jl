# PilePenetration

[![Build Status](https://github.com/KeitaNakamura/PilePenetration.jl/workflows/CI/badge.svg)](https://github.com/KeitaNakamura/PilePenetration.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/KeitaNakamura/PilePenetration.jl/branch/main/graph/badge.svg?token=DPE75E110O)](https://codecov.io/gh/KeitaNakamura/PilePenetration.jl)

## Run using [Julia](https://julialang.org) REPL

1. Add registry

```julia
pkg> add registry https://github.com/KeitaNakamura/KeitaNakamuraRegistry.git
```

2. Install package

```julia
pkg> add https://github.com/KeitaNakamura/PilePenetration.jl.git
```

3. Run simulations

```julia
julia> using PilePenetration

julia> PilePenetration.main("path/to/input.toml")
```
See files in `examples` for `input.toml`.


## Run using execution file

1. Download execution files from [releases](https://github.com/KeitaNakamura/PilePenetration.jl/releases)

2. Run the execution file

```bash
$ path/to/bin/PilePenetration[.exe] path/to/input.toml
```

See files in `examples` for `input.toml`.
If the path to `input.toml` is not given, the `input.toml` in the current directory is used.
