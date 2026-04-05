
# loading data
function data_tkachenko(i)
    # data from paper: https://doi.org/10.1016/j.jcp.2022.111901
    path_tkachenko = joinpath(@__DIR__, "data", "dam_break")
    csv_names = ["dam_break_2d_cylinder.csv", "dam_break_2d_square.csv"]

    data_path = joinpath(path_tkachenko, csv_names[i])
    data = readdlm(data_path, ','; header=false)

    x_data = data[:, 1]
    h_data = data[:, 2]

    return x_data, h_data
end

function reproduce_cylindric_dam_break_results(backend)
    # defining parameters
    reflecting_bc = Val(true)
    #backend = ROCBackend()

    tspan = (0.0, 40)
    saveat = range(tspan...; length=10)
    λ = 500
    g = 9.81

    xmin = -300.0
    xmax = 300.0
    Nx = 2728

    ymin = -300.0
    ymax = 300.0
    Ny = 2728

    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)
    dx, dy = step(gridx), step(gridy)

    ###########################################################
    # cylindric dam break 
    function initial_condtition_cylindric_dam_break(gridx, gridy, backend, reflecting_bc)
        Nx, Ny = length(gridx), length(gridy)
        dx, dy = step(gridx), step(gridy)

        # smoothing parameter
        α = 6.5

        h = zeros(Nx, Ny)
        vx = zeros(Nx, Ny)
        vy = zeros(Nx, Ny)

        # dam 
        for j in 1:Ny
            y = gridy[j]
            for i in 1:Nx
                x = gridx[i]
                r = sqrt(x^2 + y^2)
                h[i, j] = 1.0 + 0.8 * 0.5 * (1 - tanh((r - 20.0) / α))
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

    q0 = initial_condtition_cylindric_dam_break(gridx, gridy, backend, reflecting_bc)

    # Create cache
    b = zeros(Nx, Ny)
    cache = create_cache(; backend, λ, g, gridx, gridy, b, reflecting_bc)

    # Setup ODE problem
    prob = ODEProblem(rhs_split!, q0, tspan, cache)

    callback, saved_values = save_and_print_callback(saveat, save_everything=false)

    # Solve
    solve(prob, RDPK3SpFSAL35(); save_everystep=false, saveat=saveat, reltol=1e-6, abstol=1e-6, callback=callback)

    j0 = argmin(abs.(gridy)) # looking for the index of the y-value closest to 0

    begin
        h_slice = saved_values.saveval[end][:, j0]

        # data
        x_data_cylinder, h_data_cylinder = data_tkachenko(1)

        #plot
        p_cylinder = plot(gridx, h_slice, xlabel=L"x", ylabel=L"h", xlims=(40, 180), label="Hyperbolic Serre-Green-Naghdi", title="")
        scatter!(p_cylinder, x_data_cylinder, h_data_cylinder, label="numerical data (Tkachenko et al.)", right_margin=4Plots.mm)
        savefig(p_cylinder, joinpath(plots_folder, "dam_break_cylinder.pdf"))
    end

    begin

        # heatmap
        val_begin = saved_values.saveval[begin]
        val_end = saved_values.saveval[end]
        cmin = min(minimum(val_begin), minimum(val_end))
        cmax = max(maximum(val_begin), maximum(val_end))
        clims = (cmin, cmax)

        p_cylinder_heat_begin = heatmap(gridx, gridy, val_begin',
            xlabel=L"x", ylabel=L"y", title="",
            colorbar=false,
            color=:viridis,
            clims=clims,
            aspect_ratio=:equal, #colorbar_title = L"h",
            xlims=(gridx[begin], gridx[end]),
            ylims=(gridy[begin], gridy[end]))

        p_cylinder_heat_end = heatmap(gridx, gridy, val_end',
            xlabel=L"x", ylabel=L"y", title="",
            # colorbar = false, 
            color=:viridis,
            clims=clims,
            aspect_ratio=:equal, #colorbar_title = L"h",
            xlims=(gridx[begin], gridx[end]),
            ylims=(gridy[begin], gridy[end]))


        p_heat_cylinder = plot(p_cylinder_heat_begin, p_cylinder_heat_end, size=(1100, 500),
            left_margin=5Plots.mm, right_margin=0Plots.mm, top_margin=-10Plots.mm, bottom_margin=-10Plots.mm, layout=Plots.grid(1, 2, widths=[0.435, 0.565]))
        savefig(p_heat_cylinder, joinpath(plots_folder, "dam_break_cylinder_heat_2D.pdf"))
    end
end

function reproduce_square_dam_break_results(backend)
    ###########################################################
    # square dam break

    # defining parameters
    reflecting_bc = Val(true)
    #backend = ROCBackend()

    tspan = (0.0, 40)
    saveat = range(tspan...; length=10)
    λ = 500
    g = 9.81

    xmin = -300.0
    xmax = 300.0
    Nx = 2728

    ymin = -300.0
    ymax = 300.0
    Ny = 2728

    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)
    dx, dy = step(gridx), step(gridy)

    function initial_condtition_square_dam_break(gridx, gridy, backend, reflecting_bc)
        Nx, Ny = length(gridx), length(gridy)
        dx, dy = step(gridx), step(gridy)


        # Create 2D arrays
        h = zeros(Nx, Ny)
        vx = zeros(Nx, Ny)
        vy = zeros(Nx, Ny)


        # smoothing parameter
        α = 6.5
        # dam
        for j in 1:Ny
            y = gridy[j]
            for i in 1:Nx
                x = gridx[i]

                # smooth transition functions for x and y directions
                fx = 0.5 * (1 - tanh((abs(x) - 40.0) / α))
                fy = 0.5 * (1 - tanh((abs(y) - 40.0) / α))

                # in the square region defined by abs(x) < 40 and abs(y) < 40 h will be close to 1.8 outside it will be close to 1.0
                h[i, j] = 1.0 + 0.8 * fx * fy
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

    q0 = initial_condtition_square_dam_break(gridx, gridy, backend, reflecting_bc)

    # Create cache
    b = zeros(Nx, Ny)
    cache = create_cache(; backend, λ, g, gridx, gridy, b, reflecting_bc)

    # Setup ODE problem
    prob = ODEProblem(rhs_split!, q0, tspan, cache)

    callback, saved_values = save_and_print_callback(saveat, save_everything=false)
    RDPK3SpFSAL35
    # Solve
    solve(prob, RDPK3SpFSAL35(); save_everystep=false, saveat=saveat, reltol=1e-6, abstol=1e-6, callback=callback)


    j0 = argmin(abs.(gridy))

    h_slice = saved_values.saveval[end][:, j0]



    # data
    x_data_square, h_data_square = data_tkachenko(2)

    #plot
    p_square = plot(gridx, h_slice, xlabel=L"x", ylabel=L"h", label="Hyperbolic SGN", title="", size=(900, 400), legend=:top)
    scatter!(p_square, x_data_square, h_data_square, label="numerical data (Tkachenko et al.)", right_margin=5Plots.mm, alpha=0.8, markersize=5)
    savefig(p_square, joinpath(plots_folder, "dam_break_square.pdf"))

    p_square = plot(gridx, h_slice, xlabel=L"x", ylabel=L"h", xlims=(0, 220), label="Hyperbolic Serre-Green-Naghdi", title="")
    scatter!(p_square, x_data_square, h_data_square, label="numerical data (Tkachenko et al.)", right_margin=4Plots.mm)
    savefig(p_square, joinpath(plots_folder, "dam_break_square.pdf"))


    # heatmap
    p_square_heat_begin = heatmap(gridx, gridy, saved_values.saveval[begin]',
        xlabel=L"x", ylabel=L"y", title="",
        color=:cividis, aspect_ratio=:equal, colorbar_title=L"h",
        xlims=(gridx[begin], gridx[end]),
        ylims=(gridy[begin], gridy[end]))

    p_square_heat_end = heatmap(gridx, gridy, saved_values.saveval[end]',
        xlabel=L"x", ylabel=L"y", title="",
        color=:cividis, aspect_ratio=:equal, colorbar_title=L"h",
        xlims=(gridx[begin], gridx[end]),
        ylims=(gridy[begin], gridy[end]))

    p_heat_square = plot(p_square_heat_begin, p_square_heat_end, layout=(1, 2), size=(1200, 500))
    savefig(p_heat_square, joinpath(plots_folder, "dam_break_square_heat_2D.pdf"))


    # heatmap
    val_begin = saved_values.saveval[begin]
    val_end = saved_values.saveval[end]
    cmin = min(minimum(val_begin), minimum(val_end))
    cmax = max(maximum(val_begin), maximum(val_end))
    clims = (cmin, cmax)

    p_square_heat_begin = heatmap(gridx, gridy, val_begin',
        xlabel=L"x", ylabel=L"y", title="",
        colorbar=false,
        color=:viridis,
        clims=clims,
        aspect_ratio=:equal,
        xlims=(gridx[begin], gridx[end]),
        ylims=(gridy[begin], gridy[end]))

    p_square_heat_end = heatmap(gridx, gridy, val_end',
        xlabel=L"x", ylabel=L"y", title="",
        color=:viridis,
        clims=clims,
        aspect_ratio=:equal,
        xlims=(gridx[begin], gridx[end]),
        ylims=(gridy[begin], gridy[end]))

    p_heat_square = plot(p_square_heat_begin, p_square_heat_end, size=(1100, 500),
        left_margin=5Plots.mm, right_margin=0Plots.mm, top_margin=-10Plots.mm, bottom_margin=-10Plots.mm,
        layout=Plots.grid(1, 2, widths=[0.435, 0.565]))
    savefig(p_heat_square, joinpath(plots_folder, "dam_break_square_heat_2D.pdf"))
end
