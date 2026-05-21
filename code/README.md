# Reproducing Numerical Results

This document provides instructions for reproducing the numerical experiments and benchmarking results from the paper.

**Note:** The `backend` parameter can be one of the following:

- `CPU()` for CPU execution (either x86 or ARM)
- `CUDABackend()` for Nvidia GPUs (requires CUDA-capable hardware)
- `ROCBackend()` for AMD GPUs (requires ROCm-capable hardware)
- `oneAPIBackend()` for Intel GPUs (not tested)
- `MetalBackend()` is currently not supported as the code uses Float64 precision (though this could be changed with a few lines of code to support Metal)

## Approximate Execution Times

Approximate execution times (in seconds) for all numerical experiments, including analysis and plotting. Entries marked "OOM" indicates an out-of-memory failure.

| Experiment                              | CPU Intel i7-1185G7 | Apple M4    | AMD MI210 | NVIDIA H200  |
|-----------------------------------------|---------------------| ------------|-----------|--------------|
| Hyperbolic Soliton                      | 37.8                | 2.7         | 5.6       | 14.8         |
| One-dimensional Solitary Wave           | 102.8               | 28.0        | 114.3     | 44.7         |
| Manufactured Solution                   | 10057.9             | 2962.3      | 170.7     | 237.0        |
| Semidiscrete Energy Conservation        | OOM                 | OOM         | 691.9     | 384.1        |
| Dingemans Experiment                    | 1261.5              | 288.7       | 221.7     | 19.9         |
| Head-on Collision of Solitary Waves     | 13.7                | 3.6         | 4.6       | 1.6          |
| Riemann Problem                         | 19.9                | 4.8         | 8.6       | 3.2          |
| Favre Waves (Short-time Evolution)      | 16.3                | 4.6         | 13.0      | 2.8          |
| Favre Waves (Froude Number Analysis)    | 171.4               | 40.8        | 330.4     | 9.4          |
| Reflection of Solitary Waves (A=0.075)  | 13.2                | 2.5         | 8.3       | 1.8          |
| Reflection of Solitary Waves (A=0.65)   | 8.8                 | 1.1         | 7.3       | 0.6          |
| Solitary Wave over a Gaussian Obstacle  | 1363.3              | 268.7       | 30.0      | 12.7         |
| Semi-circular Shoal (case a)            | 4670.0              | 1218.1      | 162.1     | 53.9         |
| Semi-circular Shoal (case b)            | 5045.8              | 1010.3      | 118.5     | 46.3         |
| Semi-circular Shoal (case b)            | 7128.7              | 1466.8      | 167.7     | 65.6         |
| Cylindric Dam Break                     | 3416.5              | 849.6       | 89.1      | 34.9         |
| Square Dam Break                        | 3053.9              | 882.6       | 85.5      | 40.6         |
| Long-term Cylindrical Dam Break         | OOM                 | OOM         | OOM       | 2270.0       |

## Getting Started

To reproduce the results, first include the main code file and set your desired backend:

```julia
include("code.jl")
backend = CUDABackend() # or backend = CPU(), backend = ROCBackend()
```

Then you can run any of the following commands to reproduce specific experiments.

## 5 Computing Solitary Waves Numerically
```julia
hyperbolic_soliton_results(backend)
```

## 6.1 Convergence Studies

### 6.1.1 One-Dimensional Solitary Wave Solutions
```julia
convergence_study_1D_plane_wave_combined(backend)
```

### 6.1.2 Manufactured Solution for the Two-Dimensional System
```julia
convergence_study_manufactured_solution(backend)
```

## 6.2 Semi-Discrete Energy Conservation
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

## 6.12 Performance Comparison

The benchmarking plots comparing performance across different hardware (Nvidia H200, AMD MI210, Intel i7-11th gen CPU) cannot be reproduced by simply executing the code, as they require data from runs on the specific hardware platforms. The code used to generate these plots and the collected timing data can be found in `benchmarking.jl`.
