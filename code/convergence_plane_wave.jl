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

function convergence_study_1D_plane_wave(backend)
    reflecting_bc = Val(false)

    # run simultion
    Ns = [128, 256, 512, 1024]
    directions = [:x, :y]
    L2s_pw = zeros(2, length(Ns), length(directions))
    for (j, direction) in enumerate(directions)
        for (i, N) in enumerate(Ns)
            println("Running with direction in $(direction) and N = $(N)")
            saved_values, gridx, gridy, params = run_plan_wave_simulation(N, direction, backend)
            qend = setup_solitary_wave_2D(saved_values.t[end], gridx, gridy, CPU(), reflecting_bc; params...)
            for k in 1:2
                # when direction == :y, look at v_y
                # because v_x is constant 0 and vise versa
                if direction == :y && k == 2
                    k_v = k + 1
                else
                    k_v = k
                end
                L2s_pw[k, i, j] = calculate_1D_l2(saved_values.saveval[end][:, :, k_v], qend[:, :, k_v], gridx, gridy, direction)
            end
        end
    end

    # plot and save figures
    xtick = (Ns, string.(Ns))
    for dir_idx in 1:length(directions)

        if dir_idx == 1
            # title = "Convergence study for plane wave in x-direction"
            title = ""
            label2 = L"u" * " — velocity in " * L"x"
            save_as = joinpath(plots_folder, "plane_wave_convergence_study_in_x.pdf")

        else
            # title = "Convergence study for plane wave in y-direction"
            title = ""
            label2 = L"v" * " — velocity in " * L"y"
            save_as = joinpath(plots_folder, "plane_wave_convergence_study_in_y.pdf")
        end

        p_pw = plot(Ns, L2s_pw[1, :, dir_idx],
            title=title,
            label=L"h" * " — water height", marker=:o,
            xlabel=L"\mathcal{N}", ylabel=L"\Vert f_{\mathrm{num}} - f_{\mathrm{ana}} \; \Vert_{L^2}",
            yscale=:log10, xscale=:log10,
            xticks=xtick, ms=3.5,
            legend=:topright, legend_column=1,
            size=(600, 400))
        plot!(p_pw, Ns, L2s_pw[2, :, dir_idx], label=label2, marker=:o, ms=3.5)

        ylim = ylims(p_pw)
        xlim = xlims(p_pw)
        Ns_ref = range(0.8 * Ns[1], 1.2 * Ns[end], length=100)
        ref_line = Ns[1]^2.0 * mean(L2s_pw[:, 1, dir_idx]) .* (Ns_ref .^ -2.0)
        plot!(p_pw, Ns_ref, ref_line, label=L"2^{\mathrm{nd}}" * "order reference line", linestyle=:dash, color=6, xlims=xlim, ylims=ylim)

        plot!(p_pw, left_margin=2Plots.mm, bottom_margin=1Plots.mm,)

        @info "figure saved at:" savefig(p_pw, joinpath(plots_folder, save_as))

    end

end
