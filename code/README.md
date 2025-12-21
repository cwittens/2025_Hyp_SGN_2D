# Reproducing Numerical Results

This document provides instructions for reproducing the numerical experiments and benchmarking results from the paper.

**Note:** The `backend` parameter can be one of the following:

- `CPU()` for CPU execution
- `CUDABackend()` for Nvidia GPUs (requires CUDA-capable hardware)
- `ROCBackend()` for AMD GPUs (requires ROCm-capable hardware)
- `oneAPIBackend()` for Intel GPUs (not tested)
- `MetalBackend()` is currently not supported as the code uses Float64 precision (though this could be changed with a few lines of code to support Metal)

## Approximate Execution Times

The following table shows approximate execution times for each experiment (including compilation):

| Experiment | H100 (seconds) | CPU (seconds) |
|------------|----------------|---------------|
| One-dimensional Solitary Wave | 15.4 | 148.2 |
| Manufactured Solution | 230.1 | too long |
| Conservation of Energy in the Semidiscretization | 430.1 | too long |
| Dingemans Experiment | 17.9 | 2204.2 |
| Head-on Collision of Solitary Waves | 1.5 | 15.1 |
| Riemann Problem | 1.6 | 24.5 |
| Favre Waves (Short-time Evolution) | 3.0 | 38.5 |
| Favre Waves (Froude Number Analysis) | 11.9 | 259.5 |
| Reflection of Solitary Waves (A=0.075) | 1.7 | 8.6 |
| Reflection of Solitary Waves (A=0.65) | 0.8 | 5.2 |
| Solitary Wave over a Gaussian Obstacle | 12.7 | 1456.8 |

## Getting Started

To reproduce the results, first include the main code file and set your desired backend:

```julia
include("code.jl")
backend = CPU()  # or CUDABackend(), ROCBackend()
```

Then you can run any of the following commands to reproduce specific experiments.
## 5.1 Convergence Studies

### 5.1.1 One-dimensional Solitary Wave
```julia
convergence_study_1D_plane_wave(backend)
```

### 5.1.2 Manufactured Solution
```julia
convergence_study_manufactured_solution(backend)
```

## 5.2 Conservation of Energy in the Semidiscretization
```julia
reproduce_energy_variable_results(backend, (0.0, 100.0))
```
where `(0.0, 100.0)` is the time span for the simulation.

## 5.3 Dingemans Experiment
```julia
reproduce_dingemans_results(backend)
```

## 5.4 Head-on Collision of Solitary Waves
```julia
reproduce_henderson_results(backend)
```

## 5.5 Riemann Problem
```julia
reproduce_riemann_problem_results(backend)
```

## 5.6 Favre Waves

### Short-time Evolution
```julia
reproduce_favre_1D_results(backend)
```

### Froude Number Analysis
```julia
reproduce_favre_froude(backend)
```

## 5.7 Reflection of Solitary Waves from a Vertical Wall
```julia
reproduce_mitsotakis_reflecting_wave_results(backend, i)
```
where `i` is either 1 or 2 for the `A = 0.075` or `A = 0.65` case respectively.

## 5.8 Solitary Wave over a Gaussian Obstacle
```julia
reproduce_gaussian_busto_results(backend)
```

## Benchmarking Results

The benchmarking plots comparing performance across different hardware (Nvidia H100, AMD MI210, Intel i7-11th gen CPU) cannot be reproduced by simply executing the code, as they require data from runs on the specific hardware platforms. The code used to generate these plots and the collected timing data can be found in `benchmarking.jl`.


