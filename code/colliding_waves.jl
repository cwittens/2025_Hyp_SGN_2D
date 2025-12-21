function data_henderson(i)
    # data originally from Prof. Diane Henderson
    # but scrapped from paper: https://doi.org/10.1017/S0956792513000168
    path_henderson = joinpath(@__DIR__, "data", "Henderson")
    csv_names = ["1850.csv", "1860.csv", "1870.csv", "1880.csv", "1892.csv", "1900.csv", "1905.csv",
        "1910.csv", "1915.csv", "1919.csv", "1933.csv", "1950.csv", "1985.csv", "2000.csv"]


    data_path = joinpath(path_henderson, csv_names[i])
    data = readdlm(data_path, ','; header=false)

    x_data = data[:, 1]
    h_data = data[:, 2]

    return x_data, h_data
end

function reproduce_henderson_results(backend)
    reflecting_bc = Val(false)

    paper_times = [1850, 1860, 1870, 1880, 1892, 1900, 1905, 1910, 1915, 1919, 1933, 1950, 1985, 2000]
    times = (paper_times .- paper_times[1]) ./ 100 # shift and convert to seconds
    tspan = (0.0, times[end])
    λ = 500
    g = 9.81
    xmin = -10.0
    xmax = 10.0
    Nx = 4000 

    ymin = 0.0
    ymax = 1.0
    Ny = 5

    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)
    dx, dy = step(gridx), step(gridy)

    b = zeros(Nx, Ny) # Bathymetry flat
    q0_wave1 = setup_solitary_wave_2D(0.0, gridx, gridy, CPU(), reflecting_bc; coord0=0.4, A=0.01077, h∞=0.05, direction=:x, b=b)
    q0_wave2 = setup_solitary_wave_2D(0.0, gridx, gridy, CPU(), reflecting_bc; coord0=1.195, A=0.01195, h∞=0.05, direction=:x, b=b)


    h = q0_wave1[:, :, 1] .+ q0_wave2[:, :, 1] .- 0.05 # subtract h∞ once to not double count
    vx = q0_wave1[:, :, 2] .- q0_wave2[:, :, 2] # second wave going left, so negative velocity
    vy = q0_wave1[:, :, 3] .+ q0_wave2[:, :, 3]

    vx_x, _ = calculate_f_x_and_f_y(vx, dx, dy, reflecting_bc)
    w = @. -h * vx_x # most of "-h * (vx_x + vy_y) + 1.5 * (vx * b_x + vy * b_y)" is zero


    q0 = similar(q0_wave1)
    q0[:, :, 1] .= h
    q0[:, :, 2] .= vx
    q0[:, :, 3] .= vy
    q0[:, :, 4] .= h
    q0[:, :, 5] .= w
    q0 = adapt(backend, q0)

    cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)

    saveat = times
    callback, saved_values = save_and_print_callback(saveat, save_everything=false, print_every_n=100)

    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

    begin
        simulation_plots = []
        for i in 1:length(times)
            x_data, h_data = data_henderson(i)

            p = plot(gridx, 100 * saved_values.saveval[i][:, 1],
                title="t = $(round(saved_values.t[i] + paper_times[1]/100, digits=2)) s",
                xlabel="", ylabel="",
                lw=6,
                legend=false,
                xlims=(-0.25, 1.75),
                ylims=(4.9, 7.7),
            )

            scatter!(p, x_data, 100 * (h_data .+ 0.05), color=2, ms=3., markerstrokewidth=0.01, alpha=0.8)
                # label="experimental data (Henderson)", linestyle=:dot, lw=6)

            push!(simulation_plots, p)
        end
        for i in [1, 4-1, 7-1, 10-1, 13-1]
            ylabel!(simulation_plots[i], "h [cm]")
        end
        for j in [2, 3+1, 5, 6+1, 8, 9+1, 11, 12+1, 14]
            yticks!(simulation_plots[j], [5.0, 5.5, 6.0, 6.5, 7.0, 7.5], [""])
        end

        for k in 12:14
            xlabel!(simulation_plots[k], "x [m]")
        end

        legend_plot = plot(legend=true, framestyle=:none, legendfontsize=17)
        yticks!(legend_plot, [1e10], [""]) # this fixes alignment issues 
        plot!(legend_plot, [], [], label="Hyperbolic\nSerre-Green-Naghdi", color=1,)
        scatter!(legend_plot, [], [], color=2,ms=3.5, markerstrokewidth=0.01,
            label=" experimental data", )



        all_plots =  [simulation_plots..., legend_plot]
        all_plots = [simulation_plots[1:2]..., legend_plot, simulation_plots[3:end]...]
        P = plot(all_plots..., size=(1200, 1600), #suptitle="Two colliding solitary waves",
            layout=@layout([a b c; d e f; g h i; j k l; m{0.33w} w g{0.33w}]),# m{0.33w} n{0.33w} g])
            left_margin=8Plots.mm, top_margin=-2Plots.mm
        )

        @info "figure saved at:" savefig(P, joinpath(plots_folder, "colliding_waves.pdf"))
    end

end

