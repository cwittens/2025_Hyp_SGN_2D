function data_busto()
    # https://doi.org/10.1007/s10915-021-01429-8
    path_busto = joinpath(@__DIR__, "data", "Busto_et_al.csv")
    all_data = readdlm(path_busto, ','; header=false)
    x_data = all_data[:, 1]
    h_data = all_data[:, 2]
    return x_data, h_data
end

function run_wave_over_gaussian_simulation(backend, reflecting_bc, tspan, n_saves=100)

    # set up parameters
    λ = 500
    g = 9.81
    xmin, xmax = -5, 35
    ymin, ymax = -10, 10
    Nx = 1601
    Ny = 801


    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)

    b = gaussian.(gridx, gridy'; x0=0.0, y0=0.0, A=0.1, σ=1.0) # Bathymetry Gaussian
    b[b.<1e-16] .= 0.0 # set small values to zero for performance
    q0 = setup_solitary_wave_2D(0.0, gridx, gridy, backend, reflecting_bc; coord0=-3.0, A=0.0365, h∞=0.2, direction=:x, b=b)

    cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)
    cache_cpu = create_cache(backend=CPU(), λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)

    saveat = range(tspan..., n_saves)



    callback, saved_values = save_and_print_callback(saveat, save_everything=true, print_every_n=1000)

    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

    return saved_values, cache_cpu
end

function reproduce_gaussian_busto_results(backend)

    # reflecting bc to have smaller domain and no problem with wave leaving
    # the left and coming back in on the right

    # or periodic bc because this is what Busto et al. used
    reflecting_bc = Val(false)
    tspan = (0.0, 12.0)
    saved_values, cache_cpu = run_wave_over_gaussian_simulation(backend, reflecting_bc, tspan)

    idx_t = argmin(abs.(saved_values.t .- 12))
    idx_y = argmin(abs.(cache_cpu.gridy .- 0.0))

    p = plot(cache_cpu.gridx, saved_values.saveval[idx_t][:, idx_y, 1],
        label="Hyperbolic Serre-Green-Naghdi",
        xlabel="x", ylabel="h",
        # title="Wave over Gaussian Bathymetry at y = 0",
        title = "",
        legend=:bottomright,
        lw=4,
        size =(600, 400)
    )

    x_data_busto, h_data_busto = data_busto()
    plot!(p, x_data_busto, h_data_busto, label="numerical data", linestyle=:dot, lw=3)

    @info "figure saved at:" savefig(p, joinpath(plots_folder, "busto_gaussian.pdf"))


end


function reproduce_energy_variable_results(backend, tspan=(0.0, 100.0))
    @info "this will take a while because the simulation runs for a long time. When running on a CPU consider reducing the tspan."
    P = []
    for reflecting_bc in [Val(false), Val(true)]
        n_saves = 1000
        saved_values, cache_cpu = run_wave_over_gaussian_simulation(backend, reflecting_bc, tspan, n_saves)

        energy_derivatives = zeros(length(saved_values.t))

        # this calculation is done on the CPU for simplicity so it can be slow.
        for i in 1:length(saved_values.t)
            energy_derivatives[i] = compute_energy_time_derivative(saved_values.saveval[i], cache_cpu, saved_values.t[i])
        end

        if reflecting_bc == Val(true)
            title = "Split form for reflecting BCs"
            saveas = "energy_time_derivative_split_form_reflecting.pdf"
        else
            title = "Split form for periodic BCs"
            saveas = "energy_time_derivative_split_form_periodic.pdf"
        end

        idx = 1
        p = plot(saved_values.t[1:idx:end], energy_derivatives[1:idx:end],
            label=L"\langle \partial_q E, \partial_t q \rangle_M",
            xlabel=L"t", ylabel=L"\langle \partial_q E, \partial_t q \rangle_M",
            title=title,
            legend=false,)
        push!(P, p)

        @info "figure saved at:" savefig(p, joinpath(plots_folder, saveas))
    end

    yticks!(P[1], [-5e-16, -2.5e-16, 0, 2.5e-16, 5e-16], ["-5 x 10^{-16}", "-2.5 x 10^{-16}", "0", "2.5 x 10^{-16}", "5 x 10^{-16}"])
    yticks!(P[2], [-5e-16, -2.5e-16, 0, 2.5e-16, 5e-16], [""])
    ylabel!(P[2], "")
    title!(P[1], "periodic boundary conditions")
    title!(P[2], "reflecting boundary conditions")
    p = plot(P..., layout=(1, 2), size=(1200, 400), ylims=(-5e-16, 5e-16), left_margin=9Plots.mm, bottom_margin=8Plots.mm, top_margin=4Plots.mm)

    @info "figure saved at:" savefig(p, joinpath(plots_folder, "energy_time_derivative_both.pdf"))
end
