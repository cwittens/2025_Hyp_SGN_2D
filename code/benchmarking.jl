# obviously the benchmarking plot can not be reproduced by just executing the file
# as it inherently requires data from runs on different hardware.
# However, here is the code which was used to generate the benchmarking results:

#=

function problem_to_benchmark_wave_over_gaussian_simulation(Nx, backend)

    reflecting_bc = Val(false)
    tspan = (0.0, 12.0)

    # set up parameters
    λ = 500
    g = 9.81
    xmin, xmax = -5, 35
    ymin, ymax = -10, 10
    Ny = Nx ÷ 2


    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)

    b = gaussian.(gridx, gridy'; x0=0.0, y0=0.0, A=0.1, σ=1.0) # Bathymetry Gaussian
    b[b.<1e-16] .= 0.0 # set small values to zero for performance
    q0 = setup_solitary_wave_2D(0.0, gridx, gridy, backend, reflecting_bc; coord0=-3.0, A=0.0365, h∞=0.2, direction=:x, b=b)

    cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)


    saveat = range(tspan..., length=2)

    callback, _ = save_and_print_callback(saveat, save_everything=false, print_every_n=1000)

    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    t = @elapsed solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)
    nf =  sol.stats.nf
    return t, nf
end

(t, stats) = problem_to_benchmark_wave_over_gaussian_simulation(128, CPU())


NX_CPU = [32, 64, 128, 256, 512, 1024, 2048]
Times_CPU = zeros(Float64, length(NX_CPU))
N_RHS = zeros(Int, length(NX_CPU))
for (i, N) in enumerate(NX_CPU)
    println("Benchmarking wave over gaussian simulation with Nx = $N on CPU")
    t, nf = problem_to_benchmark_wave_over_gaussian_simulation(N, CPU())
    Times_CPU[i] = t
    N_RHS[i] = nf
    println("Time taken: $t seconds")
    println("Number of RHS evaluations: $nf")
end


using DelimitedFiles
data = hcat(NX_CPU, Times_CPU, N_RHS)
writedlm("benchmark_results.txt", 
         [["Nx" "Time(s)" "RHS_evals"]; data], '\t')

=#

# data collected from running the above benchmarking code on different hardware:
NX_NVIDIA = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
Times_NVIDIA = [0.381993739, 0.48364305, 0.476482433, 0.347567449, 0.50768233, 2.665141048, 15.010059529, 111.256723434, 850.87730809]
N_RHS_NVIDIA = [5774, 6049, 5624, 5694, 6949, 11290, 21092, 41756, 83539]

NX_AMD = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
Times_AMD = [1.714042799, 1.743054427, 1.494827214, 2.054493093, 2.872841173, 7.707177025, 36.791337901, 261.389745079, 1940.014781778]
N_RHS_AMD = [5864, 6089, 5609, 5804, 6934, 11310, 21072, 41756, 83539]

NX_CPU = [32, 64, 128, 256, 512, 1024, 2048]
Times_CPU = [0.228672, 1.6039269, 4.2537968, 12.0992455, 74.2211993, 404.8305322, 2677.6244409]

intel_blue = colorant"#0068B5"
amd_red = colorant"#ED1C24"
nvidia_green = colorant"#76B900"

begin
#! format: noindent
ref_line_x = range(512, 8192, length=500)
ref_line_y = (ref_line_x .^ 2.7)
ref_line_y = ref_line_y ./ maximum(ref_line_y)
ref_line_y = ref_line_y .* 5 / minimum(ref_line_y)


# Custom ticks with decimals where needed
custom_ticks = [0.4, 1.7, 10, 100, 850, 2700]
custom_tick_labels = ["0.4s", "1.7s", "10s", "100s", "850s", "2700s"]
custom_tick_labels = ["0.4", "1.7", "10", "100", "850", "2700"]
custom_x_ticks = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

p = plot(NX_CPU, Times_CPU, label="CPU i7-11th", lw=3,
    color=intel_blue, marker=:diamond, size=(600, 400), left_margin=2Plots.mm)

plot!(p, NX_AMD, Times_AMD, label="AMD MI210", lw=3,
    color=amd_red, marker=:square)

plot!(p, NX_NVIDIA, Times_NVIDIA,
    label="Nvidia H100",
    color=nvidia_green, marker=:circle,
    # xlabel="Nx", 
    xlabel=L"N_x" * " (with " * L"N_{total} = N_x \times N_x/2" * ")",
    ylabel="Wall clock time (s)",
    # title="Execution Time vs Grid Size",
    xscale=:log10, yscale=:log10,
    legend=:topleft,
    lw=3,
    yticks=(custom_ticks, custom_tick_labels),
    xticks=(custom_x_ticks, string.(custom_x_ticks))
)

ylim = ylims(p)
xlim = xlims(p)

plot!(p, ref_line_x, ref_line_y, label=L"O(N_x^{2.7})", ls=:dash, color=:black, lw=2, ylims=ylim, xlims=xlim)

@info "Saving benchmarking plot..." savefig(p, joinpath(plots_folder, "benchmarking_plot.pdf"))

end


begin
#! format: noindent
Times_RHS_NVIDIA = Times_NVIDIA ./ N_RHS_NVIDIA
Times_RHS_AMD = Times_AMD ./ N_RHS_AMD

Problem_Size = NX_NVIDIA .^ 2 .÷ 2

n = 3
println(NX_NVIDIA[end-n:end])
Problem_Size = Problem_Size[end-n:end]
Times_RHS_NVIDIA = Times_RHS_NVIDIA[end-n:end] * 1e3
Times_RHS_AMD = Times_RHS_AMD[end-n:end] * 1e3


ref_line_x = range(1e5, 1e9, length=500)
ref_line_y = (ref_line_x .^ 1)
ref_line_y = ref_line_y ./ maximum(ref_line_y)
ref_line_y = ref_line_y .* 0.05 / minimum(ref_line_y)


# Custom ticks with decimals where needed
custom_ticks = [0.25, 1, 4, 16]
custom_tick_labels = ["250μs", "1ms", "4ms", "16ms"]
custom_x_ticks = Problem_Size
custom_x_ticks_labels = ["2^{$i}" for i in Int64.(log2.(Problem_Size))]
custom_x_ticks_labels = [L"2^{19}", L"2^{21}", L"2^{23}", L"2^{25}"]

p = plot(size=(600, 400), left_margin=2Plots.mm)



plot!(p, Problem_Size, Times_RHS_AMD, label="AMD MI210", lw=3,
    color=amd_red, marker=:square)

plot!(p, Problem_Size, Times_RHS_NVIDIA,
    label="Nvidia H100",
    color=nvidia_green, marker=:circle,
    # xlabel="Nx", 
    xlabel="Number of grid points: " * L"N_{total} = N_x \times N_y",#* " (problem size ~ " * L"N_x^2/2" * ")",
    ylabel="Time per RHS evaluation",
    # title="Execution Time vs Grid Size",
    xscale=:log10, yscale=:log10,
    legend=:topleft,
    lw=3,
    yticks=(custom_ticks, custom_tick_labels),
    xticks=(Problem_Size, custom_x_ticks_labels)
)

ylim = ylims(p)
xlim = xlims(p)

plot!(p, ref_line_x, ref_line_y, label=L"O(N_{total})", ls=:dash, color=:black, lw=2, ylims=ylim, xlims=xlim)

@info "Saving benchmarking plot2..." savefig(p, joinpath(plots_folder, "benchmarking_plot_2.pdf"))

end


#= analysis of scaling exponent 
NX = NX_NVIDIA
Times = Times_NVIDIA # ./ N_RHS

NX = NX_AMD
Times = Times_AMD

NX = NX_CPU
Times = Times_CPU

# NX = NX .* NX ./ 2
n = 4 # gives last n+1 points
NX[end-n:end]

log_Nx = log10.(NX[end-n:end])
log_times = log10.(Times[end-n:end])

# Linear regression in log space
A = [log_Nx ones(length(log_Nx))]
coeffs = A \ log_times
slope = coeffs[1] =#