function run_plan_wave_simulation(N, direction, backend, reflecting_bc=Val(false))

    coord0 = 0.0
    A = 0.2
    h∞ = 1.0
    g = 9.81
    λ = 30_000

    xmin = -30.0
    xmax = 30.0

    ymin = -30.0
    ymax = 30.0

    tspan = (0.0, (xmax - xmin) / sqrt(1.2 * g)) # one period for A = 0.2 and h∞ = 1.0

    if direction == :x
        Nx = N
        Ny = 5
    elseif direction == :y
        Nx = 5
        Ny = N
    else
        error("Direction must be either :x or :y")
    end

    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)


    b = zeros(Nx, Ny) # Bathymetry flat
    q0 = setup_solitary_wave_2D(0.0, gridx, gridy, backend, reflecting_bc; coord0=coord0, A=A, h∞=h∞, direction=direction, b=b)


    cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)

    saveat = range(tspan..., 2)
    callback, saved_values = save_and_print_callback(saveat, save_everything=true, print_every_n=100_000)


    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)


    return saved_values, gridx, gridy, (; coord0, A, h∞, direction, b)
end

function convergence_study_1D_plane_wave_combined(backend)
    reflecting_bc = Val(false)

    coord0 = 0.0
    A = 0.2
    h∞ = 1.0
    g = 9.81

    xmin = -30.0
    xmax = 30.0
    ymin = -30.0
    ymax = 30.0

    tspan = (0.0, (xmax - xmin) / sqrt(1.2 * g))
    Ns = [128, 256, 512, 1024]
    directions = [:x, :y]

    # SGN soliton using λ = 30_000
    L2s_std = zeros(2, length(Ns), length(directions))
    for (j, direction) in enumerate(directions)
        for (i, N) in enumerate(Ns)
            println("Standard SGN — direction=$(direction), N=$(N)")
            saved_values, gridx, gridy, params = run_plan_wave_simulation(N, direction, backend)
            qend = setup_solitary_wave_2D(saved_values.t[end], gridx, gridy, CPU(), reflecting_bc; params...)
            for k in 1:2
                k_v = (direction == :y && k == 2) ? k + 1 : k
                L2s_std[k, i, j] = calculate_1D_l2(
                    saved_values.saveval[end][:, :, k_v], qend[:, :, k_v], gridx, gridy, direction)
            end
        end
    end

    # Hyperbolic SGN soliton using λ = 50
    λ_hyp = 50
    L2s_hyp = zeros(2, length(Ns), length(directions))
    wave = build_solitary_wave(; λ=λ_hyp, h∞=h∞, g=g, A=A)

    for (j, direction) in enumerate(directions)
        for (i, N) in enumerate(Ns)
            println("Hyperbolic SGN — direction=$(direction), N=$(N)")
            Nx, Ny = direction == :x ? (N, 5) : (5, N)
            gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
            gridy = get_grid(Ny, ymin, ymax, reflecting_bc)
            b = zeros(Nx, Ny)

            q0 = evaluate_solitary_wave_2D(
                wave, gridx, gridy, tspan[1], backend, reflecting_bc;
                coord0=coord0, direction=direction, b=b)
            cache = create_cache(; backend=backend, λ=λ_hyp, g=g,
                gridx=gridx, gridy=gridy, b=b,
                reflecting_bc=reflecting_bc)
            saveat = range(tspan..., 2)
            callback, saved_values = save_and_print_callback(
                saveat; save_everything=true, print_every_n=100_000)

            prob = ODEProblem(rhs_split!, q0, tspan, cache)
            solve(prob, RDPK3SpFSAL35();
                save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

            qend = evaluate_solitary_wave_2D(
                wave, gridx, gridy, saved_values.t[end], CPU(), reflecting_bc;
                coord0=coord0, direction=direction, b=b)

            for k in 1:2
                k_v = (direction == :y && k == 2) ? k + 1 : k
                L2s_hyp[k, i, j] = calculate_1D_l2(
                    saved_values.saveval[end][:, :, k_v], qend[:, :, k_v], gridx, gridy, direction)
            end
        end
    end

    # Combined plots
    xtick = (Ns, string.(Ns))
    Ns_ref = range(0.8 * Ns[1], 1.2 * Ns[end]; length=100)

    for dir_idx in 1:length(directions)
        if dir_idx == 1
            label2 = L"u" * " — velocity in " * L"x"
            linestyle2 = :dot
            save_as = joinpath(plots_folder, "plane_wave_convergence_study_combined_in_x.pdf")
        else
            label2 = L"v" * " — velocity in " * L"y"
            linestyle2 = :dashdotdot
            save_as = joinpath(plots_folder, "plane_wave_convergence_study_combined_in_y.pdf")
        end

        p_std = plot(Ns, L2s_std[1, :, dir_idx];
            label="", marker=:circle, linestyle=:solid, ms=4.5,
            title="classical SGN solitary wave",
            xlabel=L"\mathcal{N}",
            ylabel=L"\Vert f_{\mathrm{num}} - f_{\mathrm{ana}} \; \Vert_{L^2}",
            yscale=:log10, xscale=:log10, xticks=xtick, legend=false)
        plot!(p_std, Ns, L2s_std[2, :, dir_idx]; label="", marker=:diamond, linestyle=linestyle2, ms=4.5)
        ylim_std = ylims(p_std)
        xlim_std = xlims(p_std)
        ref_std = Ns[1]^2.0 * mean(L2s_std[:, 1, dir_idx]) .* (Ns_ref .^ -2.0)
        plot!(p_std, Ns_ref, ref_std; label="", linestyle=:solid, alpha=0.7, color=:black,
            xlims=xlim_std, ylims=ylim_std)

        p_hyp = plot(Ns, L2s_hyp[1, :, dir_idx];
            label="", marker=:circle, linestyle=:solid, ms=4.5,
            title="hyperbolic SGN solitary wave",
            xlabel=L"\mathcal{N}", ylabel="",
            yscale=:log10, xscale=:log10, xticks=xtick, legend=false)
        plot!(p_hyp, Ns, L2s_hyp[2, :, dir_idx]; label="", marker=:diamond, linestyle=linestyle2, ms=4.5)
        xlim_hyp = xlims(p_hyp)
        ref_hyp = Ns[1]^2.0 * mean(L2s_hyp[:, 1, dir_idx]) .* (Ns_ref .^ -2.0)
        plot!(p_hyp, Ns_ref, ref_hyp; label="", linestyle=:solid, alpha=0.7, color=:black,
            xlims=xlim_hyp)#, ylims=ylim_std)

        legend_plot = plot(legend=:top, framestyle=:none, legendfontsize=13, label="", legend_column=3)
        plot!(legend_plot, [], []; label=L"h" * " — water height", marker=:circle, linestyle=:solid, ms=2.5)
        plot!(legend_plot, [], []; label=label2, marker=:diamond, linestyle=linestyle2, ms=2.5)
        plot!(legend_plot, [], []; label=L"2^{\mathrm{nd}}" * " order reference line",
            linestyle=:solid, alpha=0.7, color=:black)

        p_combined = plot(p_std, p_hyp, legend_plot;
            layout=@layout([a b; g{0.04h}]), size=(1200, 480),
            left_margin=8Plots.mm, bottom_margin=5Plots.mm, right_margin=2Plots.mm)

        @info "figure saved at:" savefig(p_combined, save_as)
    end
end
