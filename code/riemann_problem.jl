function reproduce_riemann_problem_results(backend)
    # Parameters
    T = 47.434
    Nx = 4000  # x-direction grid points
    Ny = 5     # y-direction grid points

    # Grid setup
    xmin, xmax = -600.0, 600.0
    ymin, ymax = 0.0, 1.0
    tspan = (0.0, T)

    # Create grids
    reflecting_bc = Val(false)
    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)


    h_initial = zeros(Nx, Ny)

    # Apply the 1D initial condition in x-direction, constant in y-direction
    for i in 1:Nx
        for j in 1:Ny
            h_initial[i, j] = 1.0 + (1.8 - 1.0) / 2 * (1 - tanh((gridx[i]) / 2.0))
        end
    end


    # flat bathymetry
    b = zeros(Nx, Ny)
    vx_initial = zeros(Nx, Ny)
    vy_initial = zeros(Nx, Ny)

    w = zeros(Nx, Ny)  # with b, vx, and vy zero, w is zero

    λ = 500.0
    g = 9.81

    cache = create_cache(
        backend=backend,
        λ=λ,
        g=g,
        gridx=gridx,
        gridy=gridy,
        b=b,
        reflecting_bc=reflecting_bc
    )



    q0 = zeros(Nx, Ny, 5)
    q0[:, :, 1] .= h_initial  # h
    q0[:, :, 2] .= vx_initial # vx
    q0[:, :, 3] .= vy_initial # vy
    q0[:, :, 4] .= h_initial  # eta (initialized same as h)
    q0[:, :, 5] .= w       # w (initialized to 0)

    q0 = adapt(backend, q0)


    saveat = range(tspan..., 2)
    callback, saved_values = save_and_print_callback(saveat, save_everything=false)
    prob = ODEProblem(rhs_split!, q0, tspan, cache)

    # Solve
    solve(prob, RDPK3SpFSAL35(), reltol=1e-6, abstol=1e-6, save_everystep=false, callback=callback)

    # Extract solution at final time
    h_final = saved_values.saveval[end]

    # any y_idx would do
    y_idx = 3

    # get cross-section at y=0
    h_final_cross = h_final[:, y_idx]
    h_initial_cross = h_initial[:, y_idx]

    # Plot


    p = plot(gridx, h_initial_cross,
        label=L"h^{0}",
        color=2,
        linestyle=:dash;
    )


    plot!(p, gridx, h_final_cross,
        xlabel="x",
        ylabel="h",
        label=L"h^{num}",
        legend=:bottomleft,
        grid=true,
        markersize=2,
        color=1,
        size=(600, 400),
        right_margin=3Plots.mm
    )



    plot!(p, [50, 300], [1.37, 1.37],
        label=L"h^* = 1.37",
        linestyle=:dash,
        color=:black;
    )

    plot!(p, [50, 300], [1.74, 1.74],
        label=L"h^m = 1.74",
        linestyle=:dash,
        color=:black;
    )

    xlims!(p, -300, 300)


    @info "figure saved at:" savefig(p, joinpath(plots_folder, "Riemann_problem.pdf"))

end


