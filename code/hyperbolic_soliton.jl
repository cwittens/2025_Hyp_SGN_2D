function hyperbolic_soliton_results(backend)

    # choose some parameters
    N = 2^10
    xmin = -40.0
    xmax = 40.0
    h∞ = 1.0
    η∞ = 1.0
    g = 9.81
    A = 0.5

    V = 2:1:8
    λ_values = 10 .^ V

    # error arrays
    errors_h_Linf = Float64[]
    errors_η_Linf = Float64[]
    errors_u_Linf = Float64[]
    errors_w_Linf = Float64[]
    residuals = Float64[]

    for λ in λ_values
        @show λ

        Φ_h_sol, Φ_η_sol, Φ_u_sol, Φ_w_sol, Φ_h_init, Φ_η_init, Φ_u_init, Φ_w_init, x, res, converged =
            solve_for_lambda(λ, N, xmin, xmax, h∞, η∞, g, A)

        # L-inf error
        error_h_Linf = norm(Φ_h_sol - Φ_h_init, Inf)
        error_η_Linf = norm(Φ_η_sol - Φ_η_init, Inf)
        error_u_Linf = norm(Φ_u_sol - Φ_u_init, Inf)
        error_w_Linf = norm(Φ_w_sol - Φ_w_init, Inf)

        push!(errors_h_Linf, error_h_Linf)
        push!(errors_η_Linf, error_η_Linf)
        push!(errors_u_Linf, error_u_Linf)
        push!(errors_w_Linf, error_w_Linf)
        push!(residuals, res)

        @show converged

    end

    # plot
    p_errors = plot(λ_values, errors_h_Linf,
        marker=:circle, linestyle=:solid, label=L"h",
        xlabel=L"\lambda",
        ylabel=L"\Vert \Phi_{\mathrm{num}} - \Phi_{\mathrm{ana}} \Vert_{\infty}",
        yscale=:log10, xscale=:log10,
        legend=:bottomleft,
        left_margin=3Plots.mm)

    plot!(p_errors, λ_values, errors_η_Linf,
        marker=:square, linestyle=:dash, label=L"\eta")
    plot!(p_errors, λ_values, errors_u_Linf,
        marker=:diamond, linestyle=:dot, label=L"u")
    plot!(p_errors, λ_values, errors_w_Linf,
        marker=:utriangle, linestyle=:dashdot, label=L"w")

    savefig(p_errors, joinpath(plots_folder, "solitary_wave_convergence.pdf"))

    #######################################
    # comparing with analytical solution 

    # parameters
    N = 2^9
    xmin = -40.0
    xmax = 40.0
    D = fourier_derivative_operator(xmin, xmax, N)
    x = SummationByPartsOperators.grid(D)
    D2 = D^2
    h∞ = 1.0
    η∞ = 1.0
    g = 9.81
    A = 0.5
    λ = 50

    #solving 
    Φ_h_sol, Φ_η_sol, Φ_u_sol, Φ_w_sol, Φ_h_init, Φ_η_init, Φ_u_init, Φ_w_init, x, res, converged =
        solve_for_lambda(λ, N, xmin, xmax, h∞, η∞, g, A)

    D = fourier_derivative_operator(xmin, xmax, N)
    x = SummationByPartsOperators.grid(D)

    p_combined = plot(x, Φ_h_sol, label=L"\Phi_h",
        color=1,
        xlabel=L"\xi", ylabel="Solution",
        title="",
        legend=:topright)
    plot!(p_combined, x, Φ_η_sol, label=L"\Phi_\eta",
        color=2)
    plot!(p_combined, x, Φ_u_sol, label=L"\Phi_u",
        color=3)
    plot!(p_combined, x, Φ_w_sol, label=L"\Phi_w",
        color=4)

    savefig(p_combined, joinpath(plots_folder, "solitary_wave_all_variables.pdf"))

    #########################################
    # compare analytical solution with solution for different λ values

    λ_vals = [50, 200]

    # solve for all λ values
    results = Dict()
    for λ in λ_vals
        @show λ
        Φ_h_sol, Φ_η_sol, Φ_u_sol, Φ_w_sol, Φ_h_init, Φ_η_init, Φ_u_init, Φ_w_init, x, res, converged =
            solve_for_lambda(λ, N, xmin, xmax, h∞, η∞, g, A)
        results[λ] = (Φ_h_sol, Φ_η_sol, Φ_u_sol, Φ_w_sol, Φ_h_init, Φ_η_init, Φ_u_init, Φ_w_init, x)
    end

    # get initial conditions (analytical solution)
    _, _, _, _, Φ_h_init, Φ_η_init, Φ_u_init, Φ_w_init, x = results[λ_vals[1]]

    # helper function for plotting
    function plot_variable(x, ana, sols, λ_vals, ylabel_str, ana_label; show_legend=true)
        p = plot(x, ana,
            label=ana_label,
            line=(3, :solid),
            color=:grey,
            xlabel=L"\xi", ylabel=ylabel_str,
            legend=show_legend ? :topright : false,
            left_margin=5Plots.mm, bottom_margin=5Plots.mm)
        colors = [1, 2]
        styles = [:dash, :dot]
        for (λ, style, color) in zip(λ_vals, styles, colors)
            plot!(p, x, sols[λ], label=L"\lambda = %$λ", line=(3, style), color=color, alpha=0.8)
        end
        return p
    end

    # get solutions
    h_sols = Dict(λ => results[λ][1] for λ in λ_vals)
    η_sols = Dict(λ => results[λ][2] for λ in λ_vals)
    u_sols = Dict(λ => results[λ][3] for λ in λ_vals)
    w_sols = Dict(λ => results[λ][4] for λ in λ_vals)

    # plot
    p_h = plot_variable(x, Φ_h_init, h_sols, λ_vals, L"\Phi_h", L"\Phi_h^{\mathrm{ana}}"; show_legend=false)
    p_η = plot_variable(x, Φ_η_init, η_sols, λ_vals, L"\Phi_\eta", L"\Phi_\eta^{\mathrm{ana}}"; show_legend=false)
    p_u = plot_variable(x, Φ_u_init, u_sols, λ_vals, L"\Phi_u", L"\Phi_u^{\mathrm{ana}}"; show_legend=false)
    p_w = plot_variable(x, Φ_w_init, w_sols, λ_vals, L"\Phi_w", L"\Phi_w^{\mathrm{ana}}"; show_legend=false)

    legend_plot = plot(legend=:top, framestyle=:none, legendfontsize=13, legend_column=3)
    plot!(legend_plot, [], [], label=L"\Phi^{\mathrm{ana}}", color=:grey, lw=3)
    for (λ, style, color) in zip(λ_vals, [:dash, :dot], [1, 2])
        plot!(legend_plot, [], [], label=L"\lambda = %$λ", line=(3, style), color=color, alpha=0.8)
    end

    final_plot = plot(p_h, p_η, p_u, p_w, legend_plot,
        layout=@layout([a b; c d; e{0.07h}]),
        size=(700, 650), left_margin=2Plots.mm, bottom_margin=-3Plots.mm, top_margin=0Plots.mm, right_margin=1Plots.mm)
    savefig(final_plot, joinpath(plots_folder, "solution_different_lambda.pdf"))

    ##############################################################
    # testing: convergence of initial condition to analytical solution for different λ values
    t = 0.0
    xmin = -40.0
    xmax = 40.0
    ymin = -40.0
    ymax = 40.0
    Nx = 2^11
    Ny = 2^11
    reflecting_bc = Val(false)
    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)

    direction = :x

    q0_ana = setup_solitary_wave_2D(t, gridx, gridy, backend, reflecting_bc, direction=direction)

    err_vals_h = []
    err_vals_vx = []
    err_vals_eta = []
    err_vals_w = []

    λ_test_vals = [100, 500, 2500, 10000, 20000]

    for λ in λ_test_vals
        @show λ
        q0 = setup_hyperbolic_solitary_wave_2D(t, gridx, gridy, backend, reflecting_bc; λ=λ, N=2^11, direction=direction)

        h = Array(q0[:, :, 1])
        h_ana = Array(q0_ana[:, :, 1])

        vx = Array(q0[:, :, 2])
        vx_ana = Array(q0_ana[:, :, 2])
        eta = Array(q0[:, :, 4])
        eta_ana = Array(q0_ana[:, :, 4])
        w = Array(q0[:, :, 5])
        w_ana = Array(q0_ana[:, :, 5])

        push!(err_vals_h, maximum(abs.(h - h_ana)))
        push!(err_vals_vx, maximum(abs.(vx - vx_ana)))
        push!(err_vals_eta, maximum(abs.(eta - eta_ana)))
        push!(err_vals_w, maximum(abs.(w - w_ana)))
    end

    p_err = plot(λ_test_vals, err_vals_h, marker=:circle, linestyle=:solid,
        label=L"h", xlabel=L"\lambda", ylabel=L"L_\infty",
        xscale=:log10, yscale=:log10,
        legend=:topright,
        size=(700, 500))

    plot!(p_err, λ_test_vals, err_vals_vx, marker=:diamond, linestyle=:dot,
        label=L"v_x")

    plot!(p_err, λ_test_vals, err_vals_eta, marker=:square, linestyle=:dashdot,
        label=L"\eta")

    plot!(p_err, λ_test_vals, err_vals_w, marker=:utriangle, linestyle=:dashdot,
        label=L"w")

    savefig(p_err, joinpath(plots_folder, "solitary_wave_initial_condition_convergence.pdf"))

end