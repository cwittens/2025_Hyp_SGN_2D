"""
[^Dingemans1994]:
    Dingemans (1994):
    Comparison of computations with Boussinesq-like models and laboratory measurements.
    [URL: https://resolver.tudelft.nl/uuid:c2091d53-f455-48af-a84b-ac86680455e9](https://resolver.tudelft.nl/uuid:c2091d53-f455-48af-a84b-ac86680455e9)

[^Dingemans1997]:
    Dingemans (1997):
    Water Wave Propagation Over Uneven Bottoms (In 2 Parts).
    [DOI: 10.1142/1241](https://doi.org/10.1142/1241)
"""
function data_dingemans()
    path_dingemans = joinpath(@__DIR__, "data", "Dingemans.csv")
    all_data, _ = readdlm(path_dingemans, ','; header=true)
    t_values = all_data[:, 1]
    x_values = (3.04, 9.44, 20.04, 26.04, 30.44, 37.04)
    experimental_data = all_data[:, 2:end]
    return t_values, x_values, experimental_data
end


function initial_condition_dingemans(x)
    g = 9.81
    h0 = 0.8
    A = 0.02
    # omega = 2*pi/(2.02*sqrt(2))
    k = 0.8406220896381442 # precomputed result of find_zero(k -> omega^2 - g * k * tanh(k * h0), 1.0) using Roots.jl
    offset = 2.4
    if x - offset < -34.5 * pi / k || x - offset > -4.5 * pi / k
        eta_prime = 0.0
    else
        eta_prime = A * cos(k * (x - offset))
    end
    v = sqrt(g / k * tanh(k * h0)) * eta_prime / h0
    if 11.01 <= x && x < 23.04
        b = 0.6 * (x - 11.01) / (23.04 - 11.01)
    elseif 23.04 <= x && x < 27.04
        b = 0.6
    elseif 27.04 <= x && x < 33.07
        b = 0.6 * (33.07 - x) / (33.07 - 27.04)
    else
        b = 0.0
    end
    h = 0.8 + eta_prime - b
    return h, v, b
end


function run_dingemans_simulation(backend, reflecting_bc)
    g = 9.81
    λ = 500
    tspan = (0.0, 70.0)
    xmin, xmax = -138.0, 46.0
    ymin, ymax = 0.0, 1.0
    # try and choose Nx so that gauges are part of the grid
    # eg N * (xmax - xmin) / 0.04 
    Nx = 23000
    Ny = 5

    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)

    dx, dy = step(gridx), step(gridy)

    h_profile = zero(gridx)
    v_profile = zero(gridx)
    b_profile = zero(gridx)

    for (i, x) in enumerate(gridx)
        h_profile[i], v_profile[i], b_profile[i] = initial_condition_dingemans(x)
    end

    h = repeat(h_profile, 1, Ny)     # Each column is the same h profile
    vx = repeat(v_profile, 1, Ny)    # Each column is the same vx profile
    b = repeat(b_profile, 1, Ny)     # Each column is the same b profile
    vy = zeros(Nx, Ny)

    b_x, b_y = calculate_f_x_and_f_y(b, dx, dy, reflecting_bc)
    vx_x, _ = calculate_f_x_and_f_y(vx, dx, dy, reflecting_bc)
    _, vy_y = calculate_f_x_and_f_y(vy, dx, dy, reflecting_bc)

    w = @. -h * (vx_x + vy_y) + 1.5 * (vx * b_x + vy * b_y)

    q0 = zeros(Nx, Ny, 5)
    q0[:, :, 1] .= h
    q0[:, :, 2] .= vx
    q0[:, :, 3] .= vy
    q0[:, :, 4] .= h
    q0[:, :, 5] .= w
    q0 = adapt(backend, q0)

    cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)

    saveat = collect(range(tspan..., 1000))
    push!(saveat, 14.0, 28.0, 42.0, 70.0) # add the times where we want to look at speficially 
    saveat = unique(sort(saveat))
    callback, saved_values = save_and_print_callback(saveat, save_everything=false, print_every_n=2000)

    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

    return saved_values, gridx, gridy, adapt(CPU(), b)
end






function reproduce_dingemans_results(backend)
    reflecting_bc = Val(false)

    saved_values, gridx, gridy, b = run_dingemans_simulation(backend, reflecting_bc)


    # plot like in https://numericalmathematics.github.io/DispersiveShallowWater.jl/stable/dingemans/#Visualization-of-Temporal-Evolution
    # curretly not saved as a figure
    #=
    times = [14.0, 28.0, 42.0, 70.0]
    y_limits = (-0.03, 0.87)
    x_limits = (-138.0, 46.0)

    snapshot_plots = []
    for time_val in times
        idx = argmin(abs.(saved_values.t .- time_val)) # get the closest point to the time_val

        p = plot(gridx, saved_values.saveval[idx][:, 1, 1] .+ b[:, 1],
            label="Numerical", title="t = $(round(saved_values.t[idx], digits=2)) s",
            ylims=y_limits,
            xlims=x_limits,
        )

        plot!(p, gridx, b[:, 1],
            label="Bathymetry",
            color=:black,
            legend=false,
        )

        push!(snapshot_plots, p)

    end

    xlims_zoom = [(-25, 15), (0, 40), (5, 45), (-100, -60)]
    snapshot_plots_zoom = [plot(snapshot_plots[i], xlims=xlims_zoom[i], ylims=(0.75, 0.85), title="Zoomed in at t = $(times[i])") for i in 1:4]


    all_plots = [snapshot_plots..., snapshot_plots_zoom...]

    plot(all_plots...,
        size=(900, 1100),
        layout=@layout([a b; c d; f g; h i]),
    )
    =#


    # compare to experimental data at gauges

    t_values, x_values, experimental_data = data_dingemans()

    tlims = [(20, 30), (25, 35), (30, 40), (35, 45), (40, 50), (45, 55)]
    snapshot_plots_time = []
    for (j, x_val) in enumerate(x_values)

        idx = argmin(abs.(gridx .- x_val)) # get the closest point to the x_val
        gauge_over_time = [saved_values.saveval[i][idx, 1] + b[idx, 1] for i in 1:length(saved_values.t)]
        p = plot(saved_values.t, gauge_over_time, xlims=tlims[j], color=1, xlabel="t")

        plot!(p, t_values, experimental_data[:, j],
            title="Gauge at x = $(x_val)",
            xlims=tlims[j],
            ylims=(0.765, 0.865),
            label="experimental data",
            linestyle=:dot, lw=3,
            color=2,
            legend=false
        )


        push!(snapshot_plots_time, p)

    end

    ylabel!(snapshot_plots_time[1], "η")
    ylabel!(snapshot_plots_time[3], "η")
    ylabel!(snapshot_plots_time[5], "η")
    yticks!(snapshot_plots_time[2], [0.78, 0.80, 0.82, 0.84, 0.86], [""])
    yticks!(snapshot_plots_time[4], [0.78, 0.80, 0.82, 0.84, 0.86], [""])
    yticks!(snapshot_plots_time[6], [0.78, 0.80, 0.82, 0.84, 0.86], [""])


    legend_plot = plot(legend=:top, framestyle=:none, legendfontsize=15)
    plot!(legend_plot, [], [], label="Hyperbolic Serre-Green-Naghdi", color=1,)
    plot!(legend_plot, [], [], label="experimental data", linestyle=:dot, lw=3, color=2)


    all_plots2 = [snapshot_plots_time..., legend_plot]
    p_gauges = plot(all_plots2..., layout=@layout([a b; c d; e f; g{0.05h}]),
        size=(900, 900), suptitle="", right_margin=5Plots.mm, bottom_margin=3Plots.mm,)

    @info "figure saved at:" savefig(p_gauges, joinpath(plots_folder, "dingemans_gauges.pdf"))

    # plot setup

    p_inital = plot([], [], color=1, label="Initial condition")
    plot!(p_inital, [], [], color=:black, label="Bathymetry")
    vline!(p_inital, collect(x_values), label="Gauge locations", color=2, linestyle=:dash, alpha=0.5)

    plot!(p_inital, gridx, saved_values.saveval[1][:, 1, 1] .+ b[:, 1],
        label="",
        color=1,
    )

    plot!(p_inital, gridx, b[:, 1],
        label="",
        color=:black,
        legend=:bottomleft,
        # title="Initial setup",
        title="",
        xlabel="x",
        ylabel="η",
        ylims=(-0.03, 0.90),
        size = (600, 400)
    )



    @info "figure saved at:" savefig(p_inital, joinpath(plots_folder, "dingemans_initial_setup.pdf"))

end