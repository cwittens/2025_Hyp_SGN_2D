function long_term_dam_break_results(backend)

    # defining parameters
    reflecting_bc = Val(true)
    #backend = ROCBackend()

    tspan = (0.0, 300.0)
    # saveat = range(tspan...; length=61)
    saveat = [30, 50, 65, 80, 90, 110, 125, 145, 160, 180, 200, 220, 240, 260, 280, 300]
    λ = 500
    g = 9.81

    xmin = -1100.0
    xmax = 1100.0
    Nx = 10_001
    ymin = -1100.0
    ymax = 1100.0
    Ny = 10_001
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


    # coursen the saved values to save memory
    function save_func_h_reduced(u, t, integrator, course)
        return copy(adapt(CPU(), u[:, :, 1][1:course:end, 1:course:end]))
    end

    function save_and_print_callback_reduced(saveat; course=1)
        # reset counter
        step_counter = Ref(0)
        print_every_n = 100
        # Callback that increments counter and prints every 100 steps
        function print_condition(u, t, integrator)
            step_counter[] += 1
            return step_counter[] % print_every_n == 0
        end

        function print_affect!(integrator)
            println("Step $(step_counter[]), t = $(integrator.t)")
            flush(stdout)
        end

        print_cb = DiscreteCallback(print_condition, print_affect!, save_positions=(false, false))

        saved_values = SavedValues(Float64, Matrix{Float64})
        save_func_h_reduced_closed(u, t, integrator) = save_func_h_reduced(u, t, integrator, course)
        save_cb = SavingCallback(save_func_h_reduced_closed, saved_values, saveat=saveat)




        return (CallbackSet(save_cb, print_cb), saved_values)

    end

    courser = 2
    callback, saved_values = save_and_print_callback_reduced(saveat, course=courser)
    # callback, saved_values = save_and_print_callback(saveat, save_everything=false)



    # Solve
    solve(prob, RDPK3SpFSAL35(); save_everystep=false, saveat=saveat, reltol=1e-6, abstol=1e-6, callback=callback)


    gridx_course = gridx[1:courser:end]
    gridy_course = gridy[1:courser:end]

    j0 = argmin(abs.(gridy_course)) # looking for the index of the y-value closest to 0
    i0 = argmin(abs.(gridx_course)) # looking for the index of the x-value closest to 0



    begin
        h_slices = [saved_values.saveval[t][:, j0] for t in 1:length(saved_values.t)]

        # wave over time
        p_over_time = plot(legend=:bottom, legend_columns=3, ylims=(0.87, 1.22))
        for t in [30, 50, 80, 110, 145, 180, 220, 260, 300]
            i = argmin(abs.(saved_values.t .- t))
            plot!(gridx_course, h_slices[i], label="t = $(round(Int, saved_values.t[i]))", xlabel=L"x", ylabel=L"h", xlims=(0, 1100), title="")
        end

        savefig(p_over_time, joinpath(plots_folder, "dam_break_long_over_time.pdf"))
    end

    begin
        # plot the 45 degree line
        xy_45_degree = Float64[]
        for i in eachindex(gridx_course)
            s = sqrt(gridx_course[i]^2 + gridy_course[i]^2)
            if gridx_course[i] < 0
                s = -s
            end
            push!(xy_45_degree, s)
        end

        h_45_degree = [saved_values.saveval[end][i, i] for i in eachindex(gridx_course)]


        p_cylinder = plot(gridx_course, h_slices[end], xlabel=L"r", ylabel=L"h", xlims=(800, 1050), label="cross section at 0°", title="", right_margin=5Plots.mm)
        plot!(xy_45_degree, h_45_degree, label="cross section at 45°", ls=:dot)
        savefig(p_cylinder, joinpath(plots_folder, "dam_break_long_cross_sections.pdf"))
    end




    begin
        peak_amps = Float64[]
        peak_positions = Float64[]
        for t in [30, 50, 65, 80, 90, 110, 125, 145, 160, 180, 200, 220, 240, 260, 280, 300]
            i = argmin(abs.(saved_values.t .- t))
            slice = saved_values.saveval[i][:, j0]
            idx = argmax(slice)
            push!(peak_amps, slice[idx])
            push!(peak_positions, abs(gridx_course[idx]))
        end



        p_amp = plot(xlabel=L"x", ylabel=L"h_{\mathrm{max}}", title="")
        scatter!(p_amp, peak_positions, peak_amps, label="leading peak")

        r0, a0 = peak_positions[1], peak_amps[1]
        r_ref = range(peak_positions[1], peak_positions[end], length=100)
        plot!(p_amp, r_ref, (a0 - 1.0) .* sqrt(r0) ./ sqrt.(r_ref) .+ 1.0,
            label=L"\propto 1/\sqrt{r}")

        savefig(p_amp, joinpath(plots_folder, "dam_break_long_amplitude_decay.pdf"))
    end


    begin

        function soliton_profile(x, p)
            A, x0 = p
            h∞ = 1.0
            ϵ = A / h∞
            κ = sqrt(3 * ϵ / (4 * h∞^2 * (1 + ϵ)))
            return h∞ .* (1.0 .+ ϵ .* sech.(κ .* (x .- x0)) .^ 2)
        end



        fit_times = [30, 80, 220, 300]
        fit_params = Dict(
            30 => (110, 135, 0.15, 130.0),
            80 => (280, 305, 0.1, 290.0),
            # 150 => (490, 520, 0.08, 510.0),
            220 => (732, 768, 0.07, 750.0),
            300 => (990, 1025, 0.06, 1007.0),
        )
        tops = []
        bots = []
        for (k, t) in enumerate(fit_times)
            i = argmin(abs.(saved_values.t .- t))
            slice = saved_values.saveval[i][:, j0]

            lo, hi, A0, x0 = fit_params[t]
            mask = lo .< gridx_course .< hi
            x_fit = collect(gridx_course[mask])
            h_fit = slice[mask]

            fit = curve_fit(soliton_profile, x_fit, h_fit, [A0, x0])
            A_fit, x0_fit = coef(fit)
            @show A_fit round(Int, x0_fit)

            x_plot = range(x_fit[1], x_fit[end], length=300)

            yl_top = k == 1 ? L"h" : ""
            yl_bot = k == 1 ? "error of fit" : ""

            # custom xticks per panel
            xt = if t == 30
                [round(Int, x0_fit) - 7, round(Int, x0_fit), round(Int, x0_fit) + 7]
            elseif t == 80
                [round(Int, x0_fit) - 8, round(Int, x0_fit), round(Int, x0_fit) + 8]
            elseif t == 220
                [round(Int, x0_fit) - 10, round(Int, x0_fit), round(Int, x0_fit) + 10]
            elseif t == 300
                [round(Int, x0_fit) - 10, round(Int, x0_fit), round(Int, x0_fit) + 10]
            end

            xlims = if t == 30
                (x0_fit - 10, x0_fit + 10)
            elseif t == 80
                (x0_fit - 11, x0_fit + 11)
            elseif t == 220
                (x0_fit - 13, x0_fit + 13)
            elseif t == 300
                (x0_fit - 13, x0_fit + 13)
            end

            # top: no x labels, just ticks
            p_top = plot(xlabel="", ylabel=yl_top, legend=false,
                title="t = $(round(Int, saved_values.t[i]))",
                xlims=xlims,
                xticks=(xt, ["" for _ in xt]))
            plot!(p_top, x_fit, h_fit, label="numerical solution")
            plot!(p_top, x_plot, soliton_profile(x_plot, coef(fit)),
                label="soliton fit", ls=:dash)

            # bottom: x labels shown
            h_fitted = soliton_profile(x_fit, coef(fit))
            p_bot = plot(x_fit, h_fit .- h_fitted,
                xlabel=L"x", ylabel=yl_bot,
                xlims=xlims,
                title="", label="", legend=false,
                xticks=xt)

            push!(tops, p_top)
            push!(bots, p_bot)
        end

        p_soliton = plot(tops..., bots...,
            layout=Plots.grid(2, 4, heights=[0.65, 0.35]),
            size=(1200, 500), right_margin=0Plots.mm, left_margin=6Plots.mm, bottom_margin=6Plots.mm)
        savefig(p_soliton, joinpath(plots_folder, "dam_break_long_soliton_fit.pdf"))
    end
end
