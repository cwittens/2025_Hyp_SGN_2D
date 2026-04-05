# Reproducing Numerical Results

This document provides instructions for reproducing the numerical experiments and benchmarking results from the paper.

**Note:** The `backend` parameter can be one of the following:

- `CPU()` for CPU execution
- `CUDABackend()` for Nvidia GPUs (requires CUDA-capable hardware)
- `ROCBackend()` for AMD GPUs (requires ROCm-capable hardware)
- `oneAPIBackend()` for Intel GPUs (not tested)
- `MetalBackend()` is currently not supported as the code uses Float64 precision (though this could be changed with a few lines of code to support Metal)

## Approximate Execution Times

The following table shows approximate execution times (in seconds) for each experiment (including compilation):

| Experiment                              | NVIDIA H200  | AMD MI210 | CPU Intel i7-1185G7 |
|-----------------------------------------|--------------|-----------|---------------------|
| Hyperbolic Soliton                      | 82.1         | 102.1     | 166.7               |
| One-dimensional Solitary Wave           | 29.0         | 82.5      | 148.2               |
| Manufactured Solution                   | 229.2        | 170.4     | too long            |
| Semidiscrete Energy Conservation        | 481.9        | 694.2     | too long            |
| Dingemans Experiment                    | 2090         | 200.6     | 2204.2              |
| Head-on Collision of Solitary Waves     | 1.8          | 10.4      | 15.1                |
| Riemann Problem                         | 2.5          | 18.6      | 24.5                |
| Favre Waves (Short-time Evolution)      | 3.0          | 21.9      | 38.5                |
| Favre Waves (Froude Number Analysis)    | 10.9         | 266.2     | 259.5               |
| Reflection of Solitary Waves (A=0.075)  | 3.7          | 15.8      | 8.6                 |
| Reflection of Solitary Waves (A=0.65)   | 0.6          | 13.6      | 5.2                 |
| Solitary Wave over a Gaussian Obstacle  | 14.0         | 40.0      | 1456.8              |
| Semi-circular Shoal (case a)            | 50.7         | 181.9     | too long            |
| Semi-circular Shoal (case b)            | 47.2         | 108.8     | too long            |
| Semi-circular Shoal (case b)            | 63.6         | 159.5     | too long            |
| Cylindric Dam Break                     | 38.3         | 87.3      | too long            |
| Square Dam Break                        | 43.3         | 90.4      | too long            |
| Long-term Cylindrical Dam Break         | 2086.3       | too long  | too long            |

## Getting Started

To reproduce the results, first include the main code file and set your desired backend:

```julia
include("code.jl")
backend = CUDABackend() # or backend = CPU(), backend = ROCBackend()
```

Then you can run any of the following commands to reproduce specific experiments.

## 5 Relative Equilibrium Structure

### 5.1 Solitary Wave Solutions for 1D HSGN
```julia
hyperbolic_soliton_results(backend)
```

## 6.1 Convergence Studies

### 6.1.1 One-dimensional Solitary Wave
```julia
convergence_study_1D_plane_wave(backend)
```

### 6.1.2 Manufactured Solution for the Two-Dimensional System
```julia
convergence_study_manufactured_solution(backend)
```

## 6.2 Semidiscrete Energy Conservation
```julia
reproduce_energy_variable_results(backend, (0.0, 100.0))
```
where `(0.0, 100.0)` is the time span for the simulation.

## 6.3 Dingemans Experiment
```julia
reproduce_dingemans_results(backend)
```

## 6.4 Head-on Collision of Solitary Waves
```julia
reproduce_henderson_results(backend)
```

## 6.5 Riemann Problem
```julia
reproduce_riemann_problem_results(backend)
```

## 6.6 Favre Waves

### 6.6.1 Evolution of the Free-Surface Elevation
```julia
reproduce_favre_1D_results(backend)
```

### 6.6.2 Comparison with Experimental Data
```julia
reproduce_favre_froude(backend)
```

## 6.7 Reflection of Solitary Waves from a Vertical Wall
```julia
reproduce_mitsotakis_reflecting_wave_results(backend, i)
```
where `i` is either 1 or 2 for the `A = 0.075` or `A = 0.65` case respectively.

## 6.8 Solitary Wave over a Gaussian Obstacle
```julia
reproduce_gaussian_busto_results(backend)
```

## 6.9 Propagation of Periodic Waves over a Semi-circular Shoal
```julia
reproduce_semi_shoal_results(backend, i)
```
where `i` is 1, 2, or 3 for the different cases.

## 6.10 Dam Break Problems
```julia
reproduce_cylindric_dam_break_results(backend)
reproduce_square_dam_break_results(backend)
```

## 6.11 Long-term Cylindrical Dam Break
```julia
long_term_dam_break_results(backend)
```

## 6.12 Benchmarking Results

The benchmarking plots comparing performance across different hardware (Nvidia H200, AMD MI210, Intel i7-11th gen CPU) cannot be reproduced by simply executing the code, as they require data from runs on the specific hardware platforms. The code used to generate these plots and the collected timing data can be found in `benchmarking.jl`.
