# EnrichedFiniteElements

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KQEFEM.github.io/EnrichedFiniteElements.jl/dev/)
[![CI](https://github.com/KQEFEM/EnrichedFiniteElements.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/KQEFEM/EnrichedFiniteElements.jl/actions/workflows/CI.yml)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**EnrichedFiniteElements.jl** is a Julia package for implementing enriched finite element methods (EFEM), with a focus on the techniques and formulations described in the following papers:

1. [**Space–time enriched finite elements for elastodynamic wave propagation**](https://link.springer.com/article/10.1007/s00366-023-01874-z)
2. [**Space-time enriched finite elements for acoustic and elastodynamic problems**](https://www.ros.hw.ac.uk/handle/10399/4887)

This package is under active development and aims to provide a practical and efficient implementation of enriched FEM for solving partial differential equations (PDEs).

The first version will be solving the Acoustic Wave Equation (original codes were MATLAB)

---

## Installation

To install `EnrichedFiniteElements.jl`, run the following in the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/KQEFEM/EnrichedFiniteElements.jl")
```
## Note
The Docker container is only tested on Apple Silicon and WSL (Although painfully slow to buuld). 
