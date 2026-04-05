#=
Semi-circular shoal experiment (Whalin, 1971)

Reproduction of Section 9.1 of Ricchiuto & Filippini (2014)
URL: https://www.math.u-bordeaux.fr/~mricchiu/RF2014.pdf
(this is what the page numbers refer to)

Three test cases with periodic wave trains (page 37):
  (a) T=1s, A=0.0195m, h₀/λ = 0.306
  (b) T=2s, A=0.0075m, h₀/λ = 0.117
  (c) T=3s, A=0.0068m, h₀/λ = 0.074
=#

function data_semi_shoal(i)
    # data scrapped from paper: https://www.math.u-bordeaux.fr/~mricchiu/RF2014.pdf
    path = joinpath(@__DIR__, "data", "circular_shoal")

    data_path_1 = joinpath(path, "T$(i)_1.csv")
    data_path_2 = joinpath(path, "T$(i)_2.csv")
    if i != 1
        data_path_3 = joinpath(path, "T$(i)_3.csv")
    end
    data1 = readdlm(data_path_1, ','; header=false)
    data2 = readdlm(data_path_2, ','; header=false)
    if i != 1
        data3 = readdlm(data_path_3, ','; header=false)
    end
    X = []
    H = []
    push!(X, data1[:, 1])
    push!(H, data1[:, 2])
    push!(X, data2[:, 1])
    push!(H, data2[:, 2])
    if i != 1
        push!(X, data3[:, 1])
        push!(H, data3[:, 2])
    end

    return X, H
end

# Bathymetry (page 37 of RF2014)
function semicircular_shoal_bathymetry(gridx, gridy)

    # Tank width: 6.096m (page 37)
    W = 6.096

    Nx, Ny = length(gridx), length(gridy)
    b = zeros(Nx, Ny)

    for j in 1:Ny
        y = gridy[j]
        G = sqrt(max(y * (W - y), 0.0))  # semi-circular function (page 37)

        for i in 1:Nx
            x = gridx[i]

            x_start = 10.67 - G  # start of slope
            x_end = 18.29 - G    # end of slope, 1:25 gradient

            if x < x_start
                z = 0.0
            elseif x < x_end
                z = (x - x_start) / 25.0
            else
                z = 0.30480
            end

            b[i, j] = z
        end
    end

    return b
end


# Source kernel: wave generation + sponge layer
# Copied per case because KernelAbstractions closures are not straightforward.

# Case A
@kernel function source_semicircular_shoal_A!(dh, dvx, dvy, deta, dw,
    @Const(h), @Const(vx), @Const(vy), @Const(eta), @Const(w),
    @Const(gridx), @Const(gridy), t, g)
    i, j = @index(Global, NTuple)
    x = gridx[i]

    # case A
    T_wave = 1.0
    A_iwg = 0.0195
    α_iwg = 1.8
    β_iwg = 0.34

    h0 = 0.4572          # still water depth [m] (page 37)
    W = 6.096           # tank width [m] (page 37)

    x_iwg = -2.0       # generator position (page 37)

    ω = 2 * π / T_wave
    C0 = sqrt(g * h0)    # shallow water celerity

    # Gaussian width (eq. 39, page 22):
    #   d_iwg = (1/2) · α_iwg · C₀ · T / √20
    d_iwg = 0.5 * α_iwg * C0 * T_wave / sqrt(20.0)

    # Amplitude scaling (eq. 39, page 22):
    #   Γ_iwg = β_iwg · C₀ · T / h₀
    Γ_iwg = β_iwg * C0 * T_wave / h0

    # Source term: -∂_t(h_iwg) added to continuity equation (eq. 36-37, page 21)
    #   h_iwg(x,t) = f_iwg(x) · A_iwg · sin(ωt)
    #   ∂_t h_iwg  = f_iwg(x) · A_iwg · ω · cos(ωt)
    f_iwg = Γ_iwg * exp(-(x - x_iwg)^2 / d_iwg^2)
    source = -f_iwg * A_iwg * ω * cos(ω * t)

    @inbounds dh[i, j] += source

    # ================================================================
    # Sponge layer (relaxation damping)
    # ================================================================
    # Alternative to the viscous sponge in Section 5.1 (pages 18-19).
    # Damps all variables towards the still-water rest state.
    # Quadratic spatial ramp for smoothness.
    #
    # Sponge lengths: 5m at each end (page 37)
    # Domain: [-10, 36]m (page 37)
    σ_max = 5.0

    # Left sponge: x ∈ [-10, -5]
    x_left_inner = -5.0
    x_left_outer = -10.0
    # Right sponge: x ∈ [31, 36]
    x_right_inner = 31.0
    x_right_outer = 36.0

    σ = 0

    if x < x_left_inner
        if x <= x_left_outer
            σ = σ_max
        else
            frac = (x_left_inner - x) / (x_left_inner - x_left_outer)
            σ = σ_max * frac^2
        end
    elseif x > x_right_inner
        if x >= x_right_outer
            σ = σ_max
        else
            frac = (x - x_right_inner) / (x_right_outer - x_right_inner)
            σ = σ_max * frac^2
        end
    end

    if σ > zero(eltype(h))
        # Compute local still water depth for the rest state target
        y = gridy[j]
        G = sqrt(max(y * (W - y), zero(eltype(h))))

        x_start = 10.67 - G
        x_end = 18.29 - G

        if x < x_start
            z = 0
        elseif x < x_end
            z = (x - x_start) / 25
        else
            z = eltype(h)(0.30480)
        end

        h_target = h0 - z

        @inbounds dh[i, j] -= σ * (h[i, j] - h_target)
        @inbounds dvx[i, j] -= σ * vx[i, j]
        @inbounds dvy[i, j] -= σ * vy[i, j]
        @inbounds deta[i, j] -= σ * (eta[i, j] - h_target)
        @inbounds dw[i, j] -= σ * w[i, j]
    end
end

# Case B
@kernel function source_semicircular_shoal_B!(dh, dvx, dvy, deta, dw,
    @Const(h), @Const(vx), @Const(vy), @Const(eta), @Const(w),
    @Const(gridx), @Const(gridy), t, g)
    i, j = @index(Global, NTuple)
    x = gridx[i]

    # case B
    T_wave = 2.0
    A_iwg = 0.0075
    α_iwg = 1.4
    β_iwg = 0.15

    h0 = 0.4572          # still water depth [m] (page 37)
    W = 6.096           # tank width [m] (page 37)

    x_iwg = -2.0       # generator position (page 37)

    ω = 2 * π / T_wave
    C0 = sqrt(g * h0)    # shallow water celerity

    # Gaussian width (eq. 39, page 22):
    #   d_iwg = (1/2) · α_iwg · C₀ · T / √20
    d_iwg = 0.5 * α_iwg * C0 * T_wave / sqrt(20.0)

    # Amplitude scaling (eq. 39, page 22):
    #   Γ_iwg = β_iwg · C₀ · T / h₀
    Γ_iwg = β_iwg * C0 * T_wave / h0

    # Source term: -∂_t(h_iwg) added to continuity equation (eq. 36-37, page 21)
    #   h_iwg(x,t) = f_iwg(x) · A_iwg · sin(ωt)
    #   ∂_t h_iwg  = f_iwg(x) · A_iwg · ω · cos(ωt)
    f_iwg = Γ_iwg * exp(-(x - x_iwg)^2 / d_iwg^2)
    source = -f_iwg * A_iwg * ω * cos(ω * t)

    @inbounds dh[i, j] += source

    # ================================================================
    # Sponge layer (relaxation damping)
    # ================================================================
    # Alternative to the viscous sponge in Section 5.1 (pages 18-19).
    # Damps all variables towards the still-water rest state.
    # Quadratic spatial ramp for smoothness.
    #
    # Sponge lengths: 5m at each end (page 37)
    # Domain: [-10, 36]m (page 37)
    σ_max = 5.0

    # Left sponge: x ∈ [-10, -5]
    x_left_inner = -5.0
    x_left_outer = -10.0
    # Right sponge: x ∈ [31, 36]
    x_right_inner = 31.0
    x_right_outer = 36.0

    σ = 0

    if x < x_left_inner
        if x <= x_left_outer
            σ = σ_max
        else
            frac = (x_left_inner - x) / (x_left_inner - x_left_outer)
            σ = σ_max * frac^2
        end
    elseif x > x_right_inner
        if x >= x_right_outer
            σ = σ_max
        else
            frac = (x - x_right_inner) / (x_right_outer - x_right_inner)
            σ = σ_max * frac^2
        end
    end

    if σ > zero(eltype(h))
        # Compute local still water depth for the rest state target
        y = gridy[j]
        G = sqrt(max(y * (W - y), zero(eltype(h))))

        x_start = 10.67 - G
        x_end = 18.29 - G

        if x < x_start
            z = 0
        elseif x < x_end
            z = (x - x_start) / 25
        else
            z = eltype(h)(0.30480)
        end

        h_target = h0 - z

        @inbounds dh[i, j] -= σ * (h[i, j] - h_target)
        @inbounds dvx[i, j] -= σ * vx[i, j]
        @inbounds dvy[i, j] -= σ * vy[i, j]
        @inbounds deta[i, j] -= σ * (eta[i, j] - h_target)
        @inbounds dw[i, j] -= σ * w[i, j]
    end
end

# Case C
@kernel function source_semicircular_shoal_C!(dh, dvx, dvy, deta, dw,
    @Const(h), @Const(vx), @Const(vy), @Const(eta), @Const(w),
    @Const(gridx), @Const(gridy), t, g)
    i, j = @index(Global, NTuple)
    x = gridx[i]

    # case C
    T_wave = 3.0
    A_iwg = 0.0068
    α_iwg = 4.0
    β_iwg = 0.23

    h0 = 0.4572          # still water depth [m] (page 37)
    W = 6.096           # tank width [m] (page 37)

    x_iwg = -2.0       # generator position (page 37)

    ω = 2 * π / T_wave
    C0 = sqrt(g * h0)    # shallow water celerity

    # Gaussian width (eq. 39, page 22):
    #   d_iwg = (1/2) · α_iwg · C₀ · T / √20
    d_iwg = 0.5 * α_iwg * C0 * T_wave / sqrt(20.0)

    # Amplitude scaling (eq. 39, page 22):
    #   Γ_iwg = β_iwg · C₀ · T / h₀
    Γ_iwg = β_iwg * C0 * T_wave / h0

    # Source term: -∂_t(h_iwg) added to continuity equation (eq. 36-37, page 21)
    #   h_iwg(x,t) = f_iwg(x) · A_iwg · sin(ωt)
    #   ∂_t h_iwg  = f_iwg(x) · A_iwg · ω · cos(ωt)
    f_iwg = Γ_iwg * exp(-(x - x_iwg)^2 / d_iwg^2)
    source = -f_iwg * A_iwg * ω * cos(ω * t)

    @inbounds dh[i, j] += source

    # ================================================================
    # Sponge layer (relaxation damping)
    # ================================================================
    # Alternative to the viscous sponge in Section 5.1 (pages 18-19).
    # Damps all variables towards the still-water rest state.
    # Quadratic spatial ramp for smoothness.
    #
    # Sponge lengths: 5m at each end (page 37)
    # Domain: [-10, 36]m (page 37)
    σ_max = 5.0

    # Left sponge: x ∈ [-10, -5]
    x_left_inner = -5.0
    x_left_outer = -10.0
    # Right sponge: x ∈ [31, 36]
    x_right_inner = 31.0
    x_right_outer = 36.0

    σ = 0

    if x < x_left_inner
        if x <= x_left_outer
            σ = σ_max
        else
            frac = (x_left_inner - x) / (x_left_inner - x_left_outer)
            σ = σ_max * frac^2
        end
    elseif x > x_right_inner
        if x >= x_right_outer
            σ = σ_max
        else
            frac = (x - x_right_inner) / (x_right_outer - x_right_inner)
            σ = σ_max * frac^2
        end
    end

    if σ > zero(eltype(h))
        # Compute local still water depth for the rest state target
        y = gridy[j]
        G = sqrt(max(y * (W - y), zero(eltype(h))))

        x_start = 10.67 - G
        x_end = 18.29 - G

        if x < x_start
            z = 0
        elseif x < x_end
            z = (x - x_start) / 25
        else
            z = eltype(h)(0.30480)
        end

        h_target = h0 - z

        @inbounds dh[i, j] -= σ * (h[i, j] - h_target)
        @inbounds dvx[i, j] -= σ * vx[i, j]
        @inbounds dvy[i, j] -= σ * vy[i, j]
        @inbounds deta[i, j] -= σ * (eta[i, j] - h_target)
        @inbounds dw[i, j] -= σ * w[i, j]
    end
end

function reproduce_semi_shoal_results(backend, case)

    if case == 1
        source_kernel! = source_semicircular_shoal_A!(backend)
        T = 1.0
        n_transient_periods = 55
        ylims_harmonic = (0.0, 0.04)
    elseif case == 2
        source_kernel! = source_semicircular_shoal_B!(backend)
        T = 2.0
        n_transient_periods = 15
        ylims_harmonic = (0.0, 0.015)
    elseif case == 3
        source_kernel! = source_semicircular_shoal_C!(backend)
        T = 3.0
        n_transient_periods = 15
        ylims_harmonic = (0.0, 0.015)
    end

    Nx = 2000
    Ny = 265

    n_sampling_periods = 15    # sample for DFT
    samples_per_period = 100
    λ = 500

    g = 9.81
    h_still = 0.4572  # [m] still water depth (page 37)

    # Domain: [-10, 36]m × [0, 6.096]m (page 37)
    xmin, xmax = -10.0, 36.0
    ymin, ymax = 0.0, 6.096

    # Reflective BC on top/bottom channel walls
    reflecting_bc = Val(true)
    # Create grid
    gridx = get_grid(Nx, xmin, xmax, reflecting_bc)
    gridy = get_grid(Ny, ymin, ymax, reflecting_bc)

    # ---- Bathymetry ----
    b = semicircular_shoal_bathymetry(gridx, gridy)
    h0 = h_still .- b

    # ---- Initial condition: still water, no waves ----
    q0 = zeros(Nx, Ny, 5)
    q0[:, :, 1] .= h0    # h
    q0[:, :, 2] .= 0     # vx
    q0[:, :, 3] .= 0     # vy
    q0[:, :, 4] .= h0    # eta (= h in relaxed state)
    q0[:, :, 5] .= 0     # w
    q0 = adapt(backend, q0)

    cache = create_cache(backend=backend, λ=λ, g=g,
        gridx=gridx, gridy=gridy, b=b, reflecting_bc=reflecting_bc,
        source_kernel=source_kernel!)

    t_transient = n_transient_periods * T
    t_total = (n_transient_periods + n_sampling_periods) * T
    tspan = (0.0, t_total)
    saveat_sampling = range(t_transient, t_total,
        length=n_sampling_periods * samples_per_period + 1)

    callback, saved_values = save_and_print_callback(saveat_sampling, save_everything=false, print_every_n=1000)

    prob = ODEProblem(rhs_split!, q0, tspan, cache)
    solve(prob, RDPK3SpFSAL35(), save_everystep=false, reltol=1e-6, abstol=1e-6, callback=callback)

    #=
    h_b = zero(b);
    p = surface(gridx, gridy, b', colorbar=false)
    anim = @animate for i in 1:100:length(saved_values.t)
        h_b .= b .+ saved_values.saveval[i] # extract h
        t = saved_values.t[i]
        @show t
        # p1 = deepcopy(p)
        p1 = surface(gridx, gridy, b', colorbar=false)#, xlims=(5, 35))
        surface!(p1, gridx, gridy, h_b', zlims=(-0.01, 0.51), title="Still water depth at t = $(round(t, digits=1)) s", xlabel="x [m]", ylabel="y [m]")
    end

    gif(anim, fps=10)
    =#

    function harmonic_amplitude(signal, times, T, k)
        ω_k = 2π * k / T
        N = length(signal)
        re = zero(eltype(signal))
        im_ = zero(eltype(signal))
        for n in 1:N
            re += signal[n] * cos(ω_k * times[n])
            im_ += signal[n] * sin(ω_k * times[n])
        end
        return 2 * sqrt(re^2 + im_^2) / N
    end

    j_center = argmin(abs.(collect(gridy) .- ymax / 2))
    h_still_centerline = h0[:, j_center]
    n_harmonics = 3

    times = saved_values.t
    y_middle_saved_values = [saved_values.saveval[i][:, j_center] for i in 1:length(saved_values.t)]

    n_x = length(y_middle_saved_values[1])
    n_t = length(times)

    harmonics = zeros(n_x, n_harmonics)

    for ix in 1:n_x
        # Time series of free surface elevation: η(t) = h(t) - h_still
        signal = [y_middle_saved_values[it][ix] - h_still_centerline[ix] for it in 1:n_t]

        for k in 1:n_harmonics
            harmonics[ix, k] = harmonic_amplitude(signal, times, T, k)
        end
    end

    X, H = data_semi_shoal(case)
    p1 = scatter(X[1], H[1], label="", color=1, marker=:square)
    scatter!(X[2], H[2], label="", color=2, marker=:circle)
    if case != 1
        scatter!(X[3], H[3], label="", color=3, marker=:utriangle)
    end

    plot!(gridx, harmonics[:, 1], xlims=(5, 25), ylims=ylims_harmonic, label="1st harmonic",
        title="",
        legend=true, color=1,
        xlabel="x [m]", ylabel="Harmonic amplitude [m]")
    plot!(gridx, harmonics[:, 2], xlims=(5, 25), label="2nd harmonic", color=2)
    plot!(gridx, harmonics[:, 3], xlims=(5, 25), label="3rd harmonic", color=3)

    savefig(p1, joinpath(plots_folder, "semi_shoal_case_$(case).pdf"))
end
