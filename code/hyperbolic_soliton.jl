function hyperbolic_soliton_results(backend)

    # choose some parameters
    N = 2^10
    xmin = -40.0
    xmax = 40.0
    h∞ = 1.0
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

        Φ_h_sol, Φ_η_sol, Φ_u_sol, Φ_w_sol, Φ_h_init, Φ_η_init, Φ_u_init, Φ_w_init, _, res, converged =
            solve_for_lambda(λ, N, xmin, xmax, h∞, g, A)

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
        marker=:circle, linestyle=:solid, color=1, label=L"h",
        xlabel=L"\lambda",
        ylabel=L"\Vert f_{\lambda} - f_{\mathrm{SGN}} \Vert_{\infty}",
        yscale=:log10, xscale=:log10,
        legend=:bottomleft,
        left_margin=3Plots.mm)

    plot!(p_errors, λ_values, errors_u_Linf,
        marker=:diamond, linestyle=:dot, color=2, label=L"u")
    plot!(p_errors, λ_values, errors_η_Linf,
        marker=:square, linestyle=:dash, color=4, label=L"\eta")
    plot!(p_errors, λ_values, errors_w_Linf,
        marker=:utriangle, linestyle=:dashdot, color=5, label=L"w")

    ylim_p = ylims(p_errors)
    xlim_p = xlims(p_errors)
    λ_ref = 10 .^ range(log10(0.8 * λ_values[1]), log10(1.2 * λ_values[end]); length=100)
    ref_line = λ_values[1] * 0.8 * errors_u_Linf[1] ./ λ_ref
    plot!(p_errors, λ_ref, ref_line,
        linestyle=:solid, alpha=0.7, color=:black,
        label=L"1^{\mathrm{st}}" * " order reference line",
        xlims=xlim_p, ylims=ylim_p)

    savefig(p_errors, joinpath(plots_folder, "solitary_wave_convergence.pdf"))
end
