using Pkg
Pkg.activate(@__DIR__)
# Pkg.instantiate() # if everything works as expected, only run this and not "Pkg.add(...)"

# add more packages if needed here 
# Pkg.add("OrdinaryDiffEqLowStorageRK")
# Pkg.add("DiffEqCallbacks")
# Pkg.add("Plots")
# Pkg.add("KernelAbstractions")
# Pkg.add("Adapt")
# Pkg.add("CUDA")
# Pkg.add("AMDGPU")
# Pkg.add("GPUArraysCore")
# Pkg.add("LsqFit")
# Pkg.add("LaTeXStrings")
# Pkg.add("Statistics")
# Pkg.add("DelimitedFiles")
# Pkg.add("Peaks")

using OrdinaryDiffEqLowStorageRK: ODEProblem, solve, RDPK3SpFSAL35
using DiffEqCallbacks: DiscreteCallback, SavingCallback, SavedValues, CallbackSet, terminate!
using KernelAbstractions: @kernel, @index, get_backend, CPU
using Adapt: adapt
using CUDA: CUDABackend
using AMDGPU: ROCBackend
using GPUArraysCore: @allowscalar
using LsqFit: curve_fit, coef
using LaTeXStrings: @L_str
using Statistics: mean
using DelimitedFiles: readdlm
using Peaks: findmaxima
using Plots #: plot, plot!, scatter!, scatter, layout, savefig, xlabel, ylabel, title, legend, lw, colorbar!, heatmap!, ylims
include(joinpath(@__DIR__, "helper_functions.jl"))

plots_folder = joinpath(@__DIR__, "plots")
mkpath(plots_folder) # create the plots folder if it does not exist

# load all individual simulations from separate files
# backend = CPU() # or CUDABackend(), ROCBackend()

include(joinpath(@__DIR__, "convergence_plane_wave.jl"))
# convergence_study_1D_plane_wave(backend)

include(joinpath(@__DIR__, "convergence_manufactured.jl"))
# convergence_study_manufactured_solution(backend)

include(joinpath(@__DIR__, "dingemans.jl"))
# reproduce_dingemans_results(backend)

include(joinpath(@__DIR__, "colliding_waves.jl"))
# reproduce_henderson_results(backend)

include(joinpath(@__DIR__, "wave_over_gaussian_bump.jl"))
# reproduce_gaussian_busto_results(backend)
# reproduce_energy_variable_results(backend, (0.0, 100.0)) # takes a while (also because some of the analysis is being run on CPU)

include(joinpath(@__DIR__, "reflecting_wave.jl"))
# reproduce_mitsotakis_reflecting_wave_results(backend, 1)
# reproduce_mitsotakis_reflecting_wave_results(backend, 2)

include(joinpath(@__DIR__, "riemann_problem.jl"))
# reproduce_riemann_problem_results(backend)

include(joinpath(@__DIR__, "favre_waves_1D.jl"))
# reproduce_favre_1D_results(backend)

include(joinpath(@__DIR__, "favre_froude_1D.jl"))
# reproduce_favre_froude(backend)