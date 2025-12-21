function data_mitsotakis(i::Int)
    # data scrapped from numerical experiments (Figure 6) in:
    # https://doi.org/10.1002/fld.4293
    if i == 1
        times = [0.0, 45.0, 48.0, 53.0, 90.0]
    elseif i == 2
        # paper has t = 41 but this is most likely a typo
        # to here we use t = 42
        times = [0.0, 38.0, 40.0, 42.0, 70.0]
    else
        error("i must be 1 or 2")
    end

    path_folder = joinpath(@__DIR__, "data", "Mitsotakis_et_al",)
    paths_data = ["6_$(i)_dotted.csv", "6_$(i)_red.csv", "6_$(i)_green.csv", "6_$(i)_blue.csv", "6_$(i)_black.csv"]
    all_data = []
    for (i, path) in enumerate(paths_data)
        path = joinpath(path_folder, path)
        data = readdlm(path, ','; header=false)
        push!(all_data, data)
    end

    return times, all_data
end


function reproduce_mitsotakis_reflecting_wave_results(backend, k::Int)
    if k == 1
        A = 0.075
        xlim_plot = (-12, 0)
        times, data = data_mitsotakis(1)
    elseif k == 2
        A = 0.65
        xlim_plot = (-8, 0)
        times, data = data_mitsotakis(2)
    else
        error("Int must be 1 or 2")
    end

    

    reflecting_bc = Val(true)
    λ = 500
    xmin, xmax = -100.0, 0.0
    Nx = 1001
    ymin, ymax = 0.0, 1.0
    Ny = 5

    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)

    tspan = (0.0, 90.0) .* sqrt(1 / 9.81) # dimensionless time converted to seconds


    b = zeros(Nx, Ny) # Bathymetry flat
    q0 = setup_solitary_wave_2D(0.0, gridx, gridy, backend, reflecting_bc; coord0=-50.0, A=A, h∞=1.0, direction=:x, b=b)

    cache = create_cache(backend=backend, λ=λ, g=9.81, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc)

    saveat = times .* sqrt(1 / 9.81) 

    callback, saved_values = save_and_print_callback(saveat, save_everything=false, print_every_n=500)

    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

    colors = [:gray, 2, 3, 1, :black]

    begin
        p1 = scatter([], [], xlims=(-70, 0),
            xlabel="x", ylabel="h",
            color=:red,
            ms=1, alpha=0.3,
            legend_column=2, 
            label="numerical data"# from Mitsotakis et al.",
        )
        # for i in 1:length(saved_values.t)
        for i in [1, 2, 3, 4, 5]

            if i == 1
                hi = saved_values.saveval[i][:, 1]
                plot!(p1, gridx, hi, label="t* = $(times[i])", lw=5, color=:gray, ls=:dot, alpha=0.5)
            else

                hi = saved_values.saveval[i][:, 1] 
                plot!(p1, gridx, hi, label="t* = $(times[i])", lw=3, color=colors[i])
                scatter!(p1, data[i][:, 1], data[i][:, 2] .+ 1, #  +1 because paper uses b = -1 and looks at \eta
                    label="",
                    # label="t = $(times[i])", 
                    alpha=0.3, color=colors[i])
            end
        end

        current()

    end



    p2 = plot(p1, xlims=xlim_plot, ylabel="", legend=false)


    P = plot(p1, p2, layout=(1, 2), size=(1100, 400),
        bottom_margin=8Plots.mm, left_margin=8Plots.mm,
        # suptitle="Reflecting solitary wave with A = $(A)",
        suptitle="",
        )

    savefig(P, joinpath(plots_folder, "reflecting_wave_A_$(A).pdf"))
end
