begin # define source terms

    #= created using Symbolics.jl
    using Symbolics

    @variables x y t
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    @variables g p b pi λ


    b = 8 // 100 * (cos(2 * pi * x) * cos(2 * pi * y) + 1 // 2 * cos(4 * pi * x) * cos(4 * pi * y))
    h = 2 + 1 // 2 * sin(2pi * x) * sin(2pi * y) * cos(2pi * t) - b
    vx = 3 // 10 * sin(2pi * x) * sin(2pi * t)
    vy = 3 // 10 * sin(2pi * y) * sin(2pi * t)
    eta = h
    w = -h * (Dx(vx) + Dy(vy)) + 3 // 2 * (vx * Dx(b) + vy * Dy(b))


    p = λ / 3 * (eta / h) * (1 - eta / h)

    sh = Dt(h) + Dx(h * vx) + Dy(h * vy)
    h_t = -(Dx(h * vx) + Dy(h * vy))
    svx = Dt(vx) + (Dx(h * vx^2 + 1 // 2 * g * h^2 + h * p) + Dy(h * vx * vy) + (g * h + 3 // 2 * h / eta * p) * Dx(b) + h_t * vx) / h
    svy = Dt(vy) + (Dx(h * vx * vy) + Dy(h * vy^2 + 1 // 2 * g * h^2 + h * p) + (g * h + 3 // 2 * h / eta * p) * Dy(b) + h_t * vy) / h
    seta = Dt(eta) + (Dx(h * eta * vx) + Dy(h * eta * vy) + 3 // 2 * h * vx * Dx(b) + 3 // 2 * h * vy * Dy(b) - h * w + h_t * eta) / h
    sw = Dt(w) + (Dx(h * w * vx) + Dy(h * w * vy) - λ * (1 - eta / h) + h_t * w) / h
    =#

    b_initial_source(x, y) = 0.08 * (cos(2 * pi * x) * cos(2 * pi * y) + 0.5 * cos(4 * pi * x) * cos(4 * pi * y))
    h_initial_source(x, y, t) = 2 + 0.5 * sin(2 * pi * x) * sin(2 * pi * y) * cos(2 * pi * t) - b_initial_source(x, y)
    vx_initial_source(x, y, t) = 0.3 * sin(2 * pi * x) * sin(2 * pi * t)
    vy_initial_source(x, y, t) = 0.3 * sin(2 * pi * y) * sin(2 * pi * t)
    w_initial_source(x, y, t) = (3 / 2) * ((3 / 125) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t)) + (3 / 5) * pi * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * (cos(2pi * x) + cos(2pi * y)) * sin(2pi * t)

    function q_initial_source(gridx, gridy, t)
        Nx, Ny = length(gridx), length(gridy)
        q0 = zeros(Nx, Ny, 5)
        for i in 1:Nx, j in 1:Ny
            x, y = gridx[i], gridy[j]
            q0[i, j, 1] = h_initial_source(x, y, t)
            q0[i, j, 2] = vx_initial_source(x, y, t)
            q0[i, j, 3] = vy_initial_source(x, y, t)
            q0[i, j, 4] = h_initial_source(x, y, t)
            q0[i, j, 5] = w_initial_source(x, y, t)
        end
        return q0
    end


    # it would be more performant to use a single function for all sources
    # and pre calculate common terms (but performance is not a really issue atm)
    @inline function h_source(x, y, t, g)
        return -pi * sin(2pi * y) * sin(2pi * x) * sin(2pi * t) + (3 / 10) * (-(2 / 25) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) + pi * sin(2pi * y) * cos(2pi * x) * cos(2pi * t)) * sin(2pi * x) * sin(2pi * t) + (3 / 5) * pi * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * cos(2pi * y) * sin(2pi * t) + (3 / 10) * (-(2 / 25) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) + pi * cos(2pi * t) * cos(2pi * y) * sin(2pi * x)) * sin(2pi * y) * sin(2pi * t)
    end

    @inline function vx_source(x, y, t, g)
        return (50 / 1) * ((3 / 250) * pi * cos(2pi * t) * sin(2pi * x) + (9 / 2500) * pi * sin(4pi * x) * (sin(2pi * t)^2) + (1 / 50) * g * pi * sin(2pi * y) * cos(2pi * x) * cos(2pi * t) - (9 / 2500) * pi * cos(2pi * x) * sin(2pi * x) * (sin(2pi * t)^2))
    end

    @inline function vy_source(x, y, t, g)
        return (50 / 1) * ((3 / 250) * pi * sin(2pi * y) * cos(2pi * t) + (9 / 2500) * pi * sin(4pi * y) * (sin(2pi * t)^2) + (1 / 50) * g * pi * cos(2pi * t) * cos(2pi * y) * sin(2pi * x) - (9 / 2500) * pi * sin(2pi * y) * cos(2pi * y) * (sin(2pi * t)^2))
    end

    @inline function eta_source(x, y, t, g)
        return (50 / 1) * ((3 / 125) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 125) * pi * cos(2pi * y) * sin(2pi * t) - (1 / 50) * pi * sin(2pi * y) * sin(2pi * x) * sin(2pi * t) + (3 / 3125) * pi * (sin(2pi * y)^2) * cos(2pi * x) * sin(2pi * t) + (3 / 3125) * pi * sin(2pi * y) * cos(4pi * x) * sin(4pi * y) * sin(2pi * t) + (3 / 3125) * pi * cos(4pi * y) * sin(4pi * x) * sin(2pi * x) * sin(2pi * t) - (3 / 6250) * pi * cos(4pi * y) * cos(2pi * x) * cos(4pi * x) * sin(2pi * t) - (3 / 6250) * pi * cos(4pi * y) * cos(4pi * x) * cos(2pi * y) * sin(2pi * t) - (3 / 3125) * pi * (cos(2pi * x)^2) * cos(2pi * y) * sin(2pi * t) - (3 / 3125) * pi * cos(2pi * x) * (cos(2pi * y)^2) * sin(2pi * t) + (3 / 3125) * pi * cos(2pi * y) * (sin(2pi * x)^2) * sin(2pi * t) + (3 / 250) * pi * sin(2pi * y) * cos(2pi * x) * cos(2pi * t) * sin(2pi * x) * sin(2pi * t) + (3 / 250) * pi * sin(2pi * y) * cos(2pi * t) * cos(2pi * y) * sin(2pi * x) * sin(2pi * t))
    end

    @inline function w_source(x, y, t, g)
        return ((3 / 10) * ((3 / 2) * ((3 / 125) * (4(pi^2) * sin(2pi * y) * sin(2pi * x) + (8 / 1) * (pi^2) * sin(4pi * x) * sin(4pi * y)) * sin(2pi * x) * sin(2pi * t) + (3 / 125) * (-(8 / 1) * (pi^2) * cos(4pi * y) * cos(4pi * x) - 4(pi^2) * cos(2pi * x) * cos(2pi * y)) * sin(2pi * y) * sin(2pi * t) + (6 / 125) * pi * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * cos(2pi * y) * sin(2pi * t)) - (6 / 5) * (pi^2) * sin(2pi * y) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * sin(2pi * t) + ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * ((2 / 25) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) - pi * cos(2pi * t) * cos(2pi * y) * sin(2pi * x))) * sin(2pi * y) * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * sin(2pi * t) + (3 / 10) * ((3 / 2) * ((3 / 125) * (4(pi^2) * sin(2pi * y) * sin(2pi * x) + (8 / 1) * (pi^2) * sin(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(8 / 1) * (pi^2) * cos(4pi * y) * cos(4pi * x) - 4(pi^2) * cos(2pi * x) * cos(2pi * y)) * sin(2pi * x) * sin(2pi * t) + (6 / 125) * pi * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * cos(2pi * x) * sin(2pi * t)) - (6 / 5) * (pi^2) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t) + ((2 / 25) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) - pi * sin(2pi * y) * cos(2pi * x) * cos(2pi * t)) * ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t))) * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t) + (3 / 5) * pi * ((3 / 2) * ((3 / 125) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t)) + ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x))) * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * ((3 / 2) * ((3 / 125) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t)) + ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x))) * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * cos(2pi * y) * sin(2pi * t) + (3 / 10) * (-(2 / 25) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) + pi * sin(2pi * y) * cos(2pi * x) * cos(2pi * t)) * ((3 / 2) * ((3 / 125) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t)) + ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x))) * sin(2pi * x) * sin(2pi * t) + (-(3 / 10) * (-(2 / 25) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) + pi * sin(2pi * y) * cos(2pi * x) * cos(2pi * t)) * sin(2pi * x) * sin(2pi * t) - (3 / 5) * pi * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * cos(2pi * x) * sin(2pi * t) - (3 / 5) * pi * (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) * cos(2pi * y) * sin(2pi * t) - (3 / 10) * (-(2 / 25) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) + pi * cos(2pi * t) * cos(2pi * y) * sin(2pi * x)) * sin(2pi * y) * sin(2pi * t)) * ((3 / 2) * ((3 / 125) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t)) + ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x))) + (3 / 10) * ((3 / 2) * ((3 / 125) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * sin(2pi * t) + (3 / 125) * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * sin(2pi * x) * sin(2pi * t)) + ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x))) * (-(2 / 25) * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) + pi * cos(2pi * t) * cos(2pi * y) * sin(2pi * x)) * sin(2pi * y) * sin(2pi * t)) / (2 - (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) + (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) + (3 / 2) * ((6 / 125) * pi * (-2pi * sin(2pi * y) * cos(2pi * x) - (2 / 1) * pi * cos(4pi * x) * sin(4pi * y)) * sin(2pi * y) * cos(2pi * t) + (6 / 125) * pi * (-(2 / 1) * pi * cos(4pi * y) * sin(4pi * x) - 2pi * cos(2pi * y) * sin(2pi * x)) * cos(2pi * t) * sin(2pi * x)) + ((6 / 5) * (pi^2) * cos(2pi * x) * cos(2pi * t) + (6 / 5) * (pi^2) * cos(2pi * t) * cos(2pi * y)) * (-2 + (2 / 25) * ((1 / 2) * cos(4pi * y) * cos(4pi * x) + cos(2pi * x) * cos(2pi * y)) - (1 / 2) * sin(2pi * y) * cos(2pi * t) * sin(2pi * x)) + pi * ((3 / 5) * pi * cos(2pi * x) * sin(2pi * t) + (3 / 5) * pi * cos(2pi * y) * sin(2pi * t)) * sin(2pi * y) * sin(2pi * x) * sin(2pi * t)
    end
end

function run_manufactured_solution_simulation(N, reflecting_bc, backend; n_saveat=2, save_everything=true)
    Nx = N
    Ny = N
    g = 9.81
    λ = 500
    min_coordinate = -1.0
    max_coordinate = 1.0
    tspan = (0.0, 1.0)

    gridx = get_grid(Nx, min_coordinate, max_coordinate, reflecting_bc)
    gridy = get_grid(Ny, min_coordinate, max_coordinate, reflecting_bc)

    b = b_initial_source.(gridx, gridy')
    q0 = q_initial_source(gridx, gridy, 0.0)
    q0 = adapt(backend, q0)

    # define source kernel
    @kernel function apply_source_terms!(dh, dvx, dvy, deta, dw, h, vx, vy, eta, w, @Const(gridx), @Const(gridy), t, g)
        i, j = @index(Global, NTuple)
        x = eltype(dh)(gridx[i])
        y = eltype(dh)(gridy[j])
        dh[i, j] += h_source(x, y, t, g)
        dvx[i, j] += vx_source(x, y, t, g)
        dvy[i, j] += vy_source(x, y, t, g)
        deta[i, j] += eta_source(x, y, t, g)
        dw[i, j] += w_source(x, y, t, g)
    end
    source_kernel = apply_source_terms!(backend)
    cache = create_cache(backend=backend, λ=λ, g=g, gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc, source_kernel=source_kernel)



    saveat = range(tspan..., n_saveat)
    callback, saved_values = save_and_print_callback(saveat, save_everything=save_everything)


    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

    return saved_values, gridx, gridy
end

function convergence_study_manufactured_solution(backend)

    Ns_ms = [32, 64, 128, 256, 512]
    boundary_conditions = [Val(true), Val(false)]
    L2s_ms = zeros(5, length(Ns_ms), length(boundary_conditions))

    # run the simulation
    for (j, bc) in enumerate(boundary_conditions)
        for (i, N) in enumerate(Ns_ms)
            println("Running with reflecting_bc=$(bc) and N = $(N)")

            saved_values, gridx, gridy = run_manufactured_solution_simulation(N, bc, backend)

            q_num_end = saved_values.saveval[end]
            q_ana_end = q_initial_source(gridx, gridy, saved_values.t[end])

            for k in 1:5
                L2s_ms[k, i, j] = calculate_2D_l2(q_num_end[:, :, k], q_ana_end[:, :, k], gridx, gridy, bc)
            end
        end
    end

    # plot and save figures
    P = []
    xticks_ms = (Ns_ms, string.(Ns_ms))
    variables = [
        (L"h" * " — water height", :solid),
        (L"u" * " — velocity in " * L"x", :solid),
        (L"v" * " — velocity in " * L"y", :dash),
        (L"\eta" * " — auxiliary variable", :dash),
        (L"w" * " — auxiliary variable", :solid)
    ]

    for bc_idx in [2, 1]

        if bc_idx == 1
            title = "Manufactured solution convergence study\nfor reflecting boundary conditions"
            title = "reflecting boundary conditions"
            save_as = joinpath(plots_folder, "manufactured_solution_convergence_study_reflecting.pdf")
        else
            title = "Manufactured solution convergence study\nfor periodic boundary conditions"
            title = "periodic boundary conditions"
            save_as = joinpath(plots_folder, "manufactured_solution_convergence_study_periodic.pdf")
        end



        p_ms = plot(
            yscale=:log10, xscale=:log10, xticks=xticks_ms,
            title=title,
            xlabel=L"\mathcal{N}", ylabel=L"\Vert f_{\mathrm{num}} - f_{\mathrm{ana}} \; \Vert_{L^2}",
            legend=false,
            size=(700, 650)
        )

        for (i, (label, ls)) in enumerate(variables)
            plot!(p_ms, Ns_ms, L2s_ms[i, :, bc_idx], label=label, marker=:o, linestyles=ls, ms=1.5)
        end
        ylim = ylims(p_ms)
        xlim = xlims(p_ms)
        Ns_ms_ref = range(0.8 * Ns_ms[1], 1.2 * Ns_ms[end], length=100)
        ref_line = Ns_ms[1]^2.0 * mean(L2s_ms[1:3, 1, bc_idx]) .* (Ns_ms_ref .^ -2.0)
        plot!(p_ms, Ns_ms_ref, ref_line, label=L"2^{\mathrm{nd}}" * "order reference line", linestyle=:dot, ylims=ylim, xlims=xlim)

        # plot!(p_ms, left_margin=1Plots.mm, top_margin=3Plots.mm, bottom_margin=-15Plots.mm, right_margin=5Plots.mm)

        push!(P, p_ms)

        # @info "figure saved at:" savefig(p_ms, joinpath(plots_folder, save_as))

    end

    # to legend plot
    legend_plot = plot(legend=:top, framestyle=:none, legendfontsize=13, label="", legend_column=3,)
    for (i, (label, ls)) in enumerate(variables)
        plot!(legend_plot, [], [], label=label, marker=:o, linestyles=ls, ms=0.1)
    end
    plot!(legend_plot, [], [], label=L"2^{\mathrm{nd}}" * "order reference line", linestyle=:dot)

    ylabel!(P[2], "")

    push!(P, legend_plot)
    P_MS = plot(P..., layout=@layout([a b; g{0.07h}]), size=(1200, 480),
        # suptitle="Manufactured solution convergence study",
        suptitle="",
        left_margin=8Plots.mm, bottom_margin=5Plots.mm, right_margin=2Plots.mm,
    )

    @info "figure saved at:" savefig(P_MS, joinpath(plots_folder, "manufactured_solution_convergence_study_combined.pdf"))



end


#=
# for creating a animation of the manufactured solution
reflecting_bc = Val(false)
saved_values, gridx, gridy = run_manufactured_solution_simulation(512, reflecting_bc, backend, n_saveat = 100, save_everything=false);
b_plot = b_initial_source.(gridx, gridy')

h_num_end = saved_values.saveval[end]
h_ana_end = h_initial_source.(gridx, gridy', saved_values.t[end])
max_error = maximum(abs, h_num_end .- h_ana_end)

anim = @animate for i in 1:1:length(saved_values.t)
    zlims_plot = (-0.1, 2.8)
    zlims_error = (-1.2 * max_error, 1.2 * max_error)

    t = saved_values.t[i]
    @show t
    h_t = saved_values.saveval[i]
    p1 = surface(gridx, gridy, b_plot',
        xlabel="x", ylabel="y", zlabel="h",
        title="Man. solution with Ns=$(length(gridx))",
        size=(600, 500),
        color=:viridis, zlims=zlims_plot,
        colorbar=false  # This removes the colorbar
    )
    surface!(p1, gridx, gridy, (h_t .+ b_plot)')


    hs(x, y) = h_initial_source(x, y, t)
    h_ana = hs.(gridx, gridy')

    if reflecting_bc == Val(false)
        title = "Ana. solution periodic bc at t=$(round(t, digits=5))"
    else
        title = "Ana. solution reflecting bc at t=$(round(t, digits=5))"
    end

    # Create contour plot
    p2 = surface(gridx, gridy, b_plot',
        title=title,
        xlabel="x",
        ylabel="y",
        fill=true,
        size=(600, 500),
        color=:viridis, zlims=zlims_plot,
        colorbar=false  # This removes the colorbar
    )
    surface!(p2, gridx, gridy, (h_ana .+ b_plot)')

    # plot difference between numerical and analytical solution
    p3 = surface(gridx, gridy, (h_t .- h_ana)', zlims=zlims_error,
        title="h_num - h_ana",
        xlabel="x", ylabel="y",
        color=:viridis, size=(600, 500), colorbar=false
    )

    plot(p1, p2, p3, layout=(1, 3), size=(1500, 500))
end

gif(anim)

=#