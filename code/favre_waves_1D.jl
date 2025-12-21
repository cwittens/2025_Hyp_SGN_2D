function data_favre()
    # Load Favre wave reference fully nonlinear potential solutions
    # data from DOI: https://doi.org/10.1017/S0022112095002813
    data_dir = joinpath(@__DIR__, "data", "Favre")

    all_data = Dict()

    for ε in [0.1, 0.2, 0.3]
        for t in [50, 60, 70]
            filename = "dh$(ε)t$(t)_FNPF.txt"
            data_file = joinpath(data_dir, filename)
            data = readdlm(data_file, Float64)
            all_data[(ε, t)] = data
        end
    end

    return all_data
end

# set the initial condition for 1D Favre waves
function setup_favre_1d(gridx, gridy, ε, backend, reflecting_bc)
    Nx, Ny = length(gridx), length(gridy)
    dx, dy = step(gridx), step(gridy)

    # some parameters
    g = 9.81
    h0 = 0.2
    α = 5 * h0
    h_jump = ε * h0
    h1 = h0 + h_jump
    u_jump = sqrt(g * (h1 + h0) / (2 * h0 * h1)) * h_jump
    x0 = 0.0
    v0 = 0.0

    # Create 2D arrays
    h = zeros(Nx, Ny)
    vx = zeros(Nx, Ny)
    vy = zeros(Nx, Ny)

    # Set initial condition (constant in y-direction)
    for j in 1:Ny
        for i in 1:Nx
            x = gridx[i]
            h[i, j] = h0 + (h_jump / 2) * (1 - tanh((x - x0) / α))
            vx[i, j] = v0 + u_jump / 2 * (1 - tanh((x - x0) / α))
        end
    end

    # flat bathymetry
    b = zeros(Nx, Ny)
    b_x, b_y = calculate_f_x_and_f_y(b, dx, dy, reflecting_bc)

    # Calculate derivatives for initial w
    vx_x, _ = calculate_f_x_and_f_y(vx, dx, dy, reflecting_bc)
    _, vy_y = calculate_f_x_and_f_y(vy, dx, dy, reflecting_bc)

    w = @. -h * (vx_x + vy_y) + 1.5 * (vx * b_x + vy * b_y)


    q0 = zeros(Nx, Ny, 5)
    q0[:, :, 1] .= h
    q0[:, :, 2] .= vx
    q0[:, :, 3] .= vy
    q0[:, :, 4] .= h
    q0[:, :, 5] .= w

    return adapt(backend, q0)
end

##############################
# Favre waves 1D

function reproduce_favre_1D_results(backend)
    # Grid parameters
    Nx = 2000 #4000
    Ny = 5
    xmin = -50.0
    xmax = 50.0
    ymin = 0.0
    ymax = 1.0

    # Physical parameters
    g = 9.81
    λ = 500.0
    h0 = 0.2

    # time parameters
    desired_times = [50, 60, 70]
    saveat = desired_times .* sqrt(h0 / g)
    tspan = (0.0, maximum(saveat))

    # bc
    reflecting_bc = Val(false)  # use periodic bc

    # Setup grids
    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)


    # colors for plots
    colors = [1, 2, 3]

    plots = []

    for (idx, ε) in enumerate([0.1, 0.2, 0.3])
        # Initial condition
        q0 = setup_favre_1d(gridx, gridy, ε, backend, reflecting_bc)

        # Create cache
        b = zeros(Nx, Ny)
        cache = create_cache(; backend, λ, g, gridx, gridy, b, reflecting_bc)

        # Setup ODE problem
        prob = ODEProblem(rhs_split!, q0, tspan, cache)

        callback, saved_values = save_and_print_callback(saveat, save_everything=false)

        # Solve
        solve(prob, RDPK3SpFSAL35(); save_everystep=false, saveat=saveat, reltol=1e-6, abstol=1e-6, callback=callback)

        # plotting
        p = plot(xlabel=L"x/h_0", ylabel=L"(h-h_0)/h_0",
            title=L"\varepsilon = %$ε",
            legend=false, grid=true; default, legendfontsize=17, titlefontsize=20)

        # any y_idx would do
        y_idx = 3

        for r in 1:length(saved_values.t)
            current_time = saved_values.t[r]
            v = desired_times[r]

            # extract h at y=0 
            h = saved_values.saveval[r][:, y_idx]
            h_new = (h .- h0) ./ h0
            x_h0 = gridx ./ h0

            plot!(p, x_h0, h_new;
                label=L"\tilde{t} = %$(round(v, digits=2))",
                linewidth=4,
                color=colors[r])

            # Load reference data
            all_data_favre = data_favre()
            data = all_data_favre[(ε, v)]
            scatter!(data[:, 1], data[:, 2];
                label="", linestyle=:dash, color=colors[r], alpha= 0.6, ms = 2.5)
        end

        xlims!(45, 90)
        ylims_current = Plots.ylims(p)
        ylims!(p, (-0.01, ylims_current[2]))

        push!(plots, p)
    end

    legend_plot = scatter([], [], legend=:top, framestyle=:none, legendfontsize=13, label="numerical data   ", legend_column=4,  color=:red, ms=2, alpha=0.5)
    for (i, t) in enumerate(desired_times)
        plot!(legend_plot, [], [], label=L"\tilde{t} = %$t", linewidth=2.5, color=colors[i])
    end
    push!(plots, legend_plot)

    all_plots = plot(plots..., layout=@layout([a b c; g{0.05h}]), size=(1200, 420), left_margin=8Plots.mm, top_margin=3Plots.mm, bottom_margin=2Plots.mm, right_margin=2Plots.mm)


    @info "figure saved at:" savefig(all_plots, joinpath(plots_folder, "Favre_waves.pdf"))

end

