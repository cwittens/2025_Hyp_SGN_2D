function data_froude()
    # Load Froude number experimental data
    # Favre data from: https://api.semanticscholar.org/CorpusID:126909361
    # Treske data from: doi: 10.1080/00221689409498738
    data_dir = joinpath(@__DIR__, "data", "Froude")

    data_favre_100 = readdlm(joinpath(data_dir, "Favre_amplmax_100.txt"))
    data_favre_200 = readdlm(joinpath(data_dir, "Favre_amplmax_200.txt"))
    data_treske_80 = readdlm(joinpath(data_dir, "Treske_amplmax_80.txt"))
    data_treske_160 = readdlm(joinpath(data_dir, "Treske_amplmax_160.txt"))

    return data_favre_100, data_favre_200, data_treske_80, data_treske_160
end

# Initial condition

function setup_favre_1d(g, h0, gridx, gridy, ε, backend, reflecting_bc)
    Nx, Ny = length(gridx), length(gridy)
    dx, dy = step(gridx), step(gridy)
    
    # some parameters
    α = 5 * h0
    h_jump = ε * h0
    h1 = h0 + h_jump
    u_jump = sqrt(g * (h1 + h0)/(2*h0*h1)) * h_jump
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
            h[i, j] = h0 + (h_jump/2) * (1 - tanh((x - x0)/α))
            vx[i, j] = v0 + u_jump/2 * (1 - tanh((x - x0)/α))
        end
    end
    
    # flat bathymetry
    b = zeros(Nx, Ny)
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
    
    return adapt(backend, q0)
end

function reproduce_favre_froude(backend)
    Nx = 4000
    Ny = 5
    xmin = -150.0
    xmax = 150.0
    ymin = 0.0
    ymax = 1.0
    reflecting_bc = Val(false)  # periodic bc
    
    g = 9.81
    h0 = 0.2
    λ = 500.0
    tspan = (0.0, 350.0) .* sqrt(h0 / g) # non-dimensional time
    
    froude = Float64[]
    a_max = Float64[]
    
    
    for ε in 0.02:0.04:0.3
        # Create grids
        gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
        gridy = get_grid(Ny, ymin, ymax, reflecting_bc)
        
        # flat bathymetry
        b = zeros(Nx, Ny)
        
        # initial condition
        q0 = setup_favre_1d(g, h0, gridx, gridy, ε, backend, reflecting_bc)
        
        # Create cache
        cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, 
                            b=b, reflecting_bc=reflecting_bc)
        
        # Setup ODE problem
        prob = ODEProblem(rhs_split!, q0, tspan, cache)
        
        h_jump = ε * h0
        
        # Callback
        idx_x = findfirst(x -> x > 63.5, gridx)
        idx_y = div(Ny, 2) + 1  # middle y-position
        threshold = h0 + 1.1 * h_jump
        
        terminate_callback = DiscreteCallback(
            function (u, t, integrator)
                return @allowscalar u[idx_x, idx_y, 1] > threshold
            end,
            terminate!)
        
        Δx = step(gridx)
        @info "running next" ε Δx
        
        sol = solve(prob, RDPK3SpFSAL35();
                save_everystep = false,
                reltol=1e-6, abstol=1e-6,
                isoutofdomain = (u, p, t) -> any(isinf, u),
                callback = terminate_callback)
        
        u = adapt(CPU(), sol.u[end])
        h = u[:, :, 1]
        
        # any y_idx would do
        y_idx = 3
        h_cross = h[:, y_idx]
        
        idx = findall(x -> 10 <= x <= 70, gridx)
        h_view = view(h_cross, idx)
        h_max = maximum(h_view)
        
        push!(froude, sqrt((1 + ε) * (1 + ε / 2)))
        push!(a_max, (h_max - h0) / h0)
    end
    
    fig = plot(
        xguide = L"Froude number $\mathrm{Fr}$",
        yguide = L"Max. amplitude $a_{\mathrm{max}} / h_0$",
        legend = :topleft,
        size = (600, 400),
    )

    # Plot
    plot!(fig, froude, a_max,
        label = "numerical",
        color = 1,
        linewidth = 3);

    # Plot the experimental data
    data_favre_100, data_favre_200, data_treske_80, data_treske_160 = data_froude()
    scatter!(fig, data_favre_100[:,1], data_favre_100[:,2], label = "Favre", color = :black, marker= :utriangle, markersize = 7, linewidth = 2);
    scatter!(fig, data_favre_200[:,1], data_favre_200[:,2], label = "", color = :black, marker= :utriangle, markersize = 7, linewidth = 2);
    scatter!(fig, data_treske_80[:,1], data_treske_80[:,2], label = "Treske", color = :orange, marker= :dtriangle, markersize = 7, linewidth = 2);
    scatter!(fig, data_treske_160[:,1], data_treske_160[:,2], label = "", color = :orange, marker= :dtriangle, markersize = 7, linewidth = 2)
    
    @info "figure saved at:" savefig(fig, joinpath(plots_folder, "Froud_numbers.pdf"))

end