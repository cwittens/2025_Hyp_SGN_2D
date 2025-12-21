default(
    grid=true,
    box=:on,
    size=(700, 500),
    dpi=100,
    titlefont=font(16),
    linewidth=3, gridlinewidth=2,
    markersize=4, markerstrokewidth=0.1,
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=16, yguidefontsize=16,
    ztickfontsize=14, zguidefontsize=16,
    legendfontsize=14
)

function rhs_split!(dq, q, cache, t)

    # creates @view's
    h, vx, vy, eta, w = eachslice(q; dims=3)
    dh, dvx, dvy, deta, dw = eachslice(dq; dims=3)

    (; λ, g, Nx, Ny, dx_inv, dy_inv, b, b_x, b_y, source_kernel, reflecting_bc) = cache

    backend = get_backend(q)

    kernel = kernel_rhs_split!(backend)
    kernel(dh, dvx, dvy, deta, dw, h, vx, vy, eta, w, b, b_x, b_y, λ, g, Nx, Ny, dx_inv, dy_inv, reflecting_bc, ndrange=(Nx, Ny))


    if !(source_kernel isa Nothing)
        (; gridx, gridy) = cache
        source_kernel(dh, dvx, dvy, deta, dw, h, vx, vy, eta, w, gridx, gridy, t, g, ndrange=(Nx, Ny))
    end

    return nothing
end


@kernel function kernel_rhs_split!(dh, dvx, dvy, deta, dw, @Const(h), @Const(vx), @Const(vy), @Const(eta), @Const(w), @Const(b), @Const(b_x), @Const(b_y), λ, g, Nx, Ny, dx_inv, dy_inv, reflecting_bc::Val{reflecting_bool}) where reflecting_bool

    i, j = @index(Global, NTuple)

    half = eltype(h)(0.5)

    if reflecting_bool # reflecting boundary conditions
        # Redefine indices based on boundary conditions
        idx_left = ifelse(i == 1, 1, i - 1) # use backward difference at left boundary
        idx_right = ifelse(i == Nx, Nx, i + 1) # use forward difference at right boundary
        jdx_down = ifelse(j == 1, 1, j - 1) # use backward difference at bottom boundary
        jdx_up = ifelse(j == Ny, Ny, j + 1) # use forward difference at top boundary

        # Define scaling factors for finite differences
        dx_scale_inv = ifelse(i == 1 || i == Nx, dx_inv, half * dx_inv) # 1 / dx or 1 / (2*dx)
        dy_scale_inv = ifelse(j == 1 || j == Ny, dy_inv, half * dy_inv) # 1 / dy or 1 / (2*dy)

        # Define factors for applying boundary conditions
        boundary_factor_x = ifelse(i == 1, -1, ifelse(i == Nx, 1, 0)) # most of the time 0 and only at boundaries -1 or 1
        boundary_factor_y = ifelse(j == 1, -1, ifelse(j == Ny, 1, 0)) # most of the time 0 and only at boundaries -1 or 1

    else # periodic boundary conditions

        # always use central finite differences 
        idx_left = mod1(i - 1, Nx)  # i minus 1
        idx_right = mod1(i + 1, Nx)  # i plus 1
        jdx_down = mod1(j - 1, Ny)  # j minus 1
        jdx_up = mod1(j + 1, Ny)  # j plus 1


        dx_scale_inv = half * dx_inv
        dy_scale_inv = half * dy_inv


        boundary_factor_x = 0
        boundary_factor_y = 0
    end


    # h_t + h_x*u + h*u_x + h_y*v + h*v_y = 0 
    @inbounds dh[i, j] = (-((h[idx_right, j] - h[idx_left, j]) * dx_scale_inv * vx[i, j] +
                            h[i, j] * (vx[idx_right, j] - vx[idx_left, j]) * dx_scale_inv +
                            (h[i, jdx_up] - h[i, jdx_down]) * dy_scale_inv * vy[i, j] +
                            h[i, j] * (vy[i, jdx_up] - vy[i, jdx_down]) * dy_scale_inv) +
                          # Apply reflecting boundary conditions (constant 0 for periodic boundaries)
                          boundary_factor_x * h[i, j] * vx[i, j] * 2 * dx_scale_inv +
                          boundary_factor_y * h[i, j] * vy[i, j] * 2 * dy_scale_inv
    )


    # a slower but better readable version of "Apply reflecting boundary conditions"
    # if reflecting_bool
    #     if i == 1
    #         dh[i, j] -= h[i, j] * vx[i, j] * 2 * dx_scale_inv # h * v * 1 / left_boundary_weight(D1)
    #     elseif i == Nx
    #         dh[i, j] += h[i, j] * vx[i, j] * 2 * dx_scale_inv # h * v * 1 / right_boundary_weight(D1)
    #     end

    #     if j == 1
    #         dh[i, j] -= h[i, j] * vy[i, j] * 2 * dy_scale_inv # h * v * 1 / bottom_boundary_weight(D2)
    #     elseif j == Ny
    #         dh[i, j] += h[i, j] * vy[i, j] * 2 * dy_scale_inv # h * v * 1 / top_boundary_weight(D2)
    #     end
    # end

    @inbounds minus_inv_h_ij = -1 / h[i, j]

    # u
    @inbounds dvx[i, j] = minus_inv_h_ij * (
        # Gravitation: g(h(h+b))_x - g(h+b)h_x
        g * ((h[idx_right, j] * (h[idx_right, j] + b[idx_right, j]) - h[idx_left, j] * (h[idx_left, j] + b[idx_left, j])) * dx_scale_inv -
             (h[i, j] + b[i, j]) * (h[idx_right, j] - h[idx_left, j]) * dx_scale_inv)

        #  (1/2)h*(u²)_x - (1/2)h_x*u² + (1/2)(hu)_x*u - (1/2)h*u*u_x
        +
        half * h[i, j] * (vx[idx_right, j]^2 - vx[idx_left, j]^2) * dx_scale_inv
        -
        half * (h[idx_right, j] - h[idx_left, j]) * dx_scale_inv * vx[i, j]^2
        +
        half * (h[idx_right, j] * vx[idx_right, j] - h[idx_left, j] * vx[idx_left, j]) * dx_scale_inv * vx[i, j]
        -
        half * h[i, j] * vx[i, j] * (vx[idx_right, j] - vx[idx_left, j]) * dx_scale_inv

        #  (1/2)(huv)_y + (1/2)(-h_y*u*v + h*u_y*v - h*u*v_y)
        + half * (h[i, jdx_up] * vx[i, jdx_up] * vy[i, jdx_up] - h[i, jdx_down] * vx[i, jdx_down] * vy[i, jdx_down]) * dy_scale_inv  # (1/2)(huv)_y
        + half * (-(h[i, jdx_up] - h[i, jdx_down]) * dy_scale_inv * vx[i, j] * vy[i, j]                      # (1/2)(-h_y*u*v)
                  +
                  h[i, j] * (vx[i, jdx_up] - vx[i, jdx_down]) * dy_scale_inv * vy[i, j]                    # (1/2)(h*u_y*v)
                  -
                  h[i, j] * vx[i, j] * (vy[i, jdx_up] - vy[i, jdx_down]) * dy_scale_inv)                   # (1/2)(-h*u*v_y)

        # Dispersive Terme 
        + (λ / 6) * (eta[i, j]^2 / h[i, j]^2) * (h[idx_right, j] - h[idx_left, j]) * dx_scale_inv
        + (λ / 3) * (eta[idx_right, j] - eta[idx_left, j]) * dx_scale_inv
        -
        (λ / 3) * (eta[i, j] / h[i, j]) * (eta[idx_right, j] - eta[idx_left, j]) * dx_scale_inv
        -
        (λ / 6) * ((eta[idx_right, j]^2 / h[idx_right, j]) - (eta[idx_left, j]^2 / h[idx_left, j])) * dx_scale_inv

        # Bathymetrie (x-Richtung)
        +
        (λ / 2) * (1 - eta[i, j] / h[i, j]) * b_x[i, j]
    )

    #  v
    @inbounds dvy[i, j] = minus_inv_h_ij * (
        # Gravitation: g(h(h+b))_y - g(h+b)h_y 
        g * ((h[i, jdx_up] * (h[i, jdx_up] + b[i, jdx_up]) - h[i, jdx_down] * (h[i, jdx_down] + b[i, jdx_down])) * dy_scale_inv -
             (h[i, j] + b[i, j]) * (h[i, jdx_up] - h[i, jdx_down]) * dy_scale_inv)

        # (1/2)h*(v²)_y - (1/2)h_y*v² + (1/2)(hv)_y*v - (1/2)h*v*v_y 
        +
        half * h[i, j] * (vy[i, jdx_up]^2 - vy[i, jdx_down]^2) * dy_scale_inv
        -
        half * (h[i, jdx_up] - h[i, jdx_down]) * dy_scale_inv * vy[i, j]^2
        +
        half * (h[i, jdx_up] * vy[i, jdx_up] - h[i, jdx_down] * vy[i, jdx_down]) * dy_scale_inv * vy[i, j]
        -
        half * h[i, j] * vy[i, j] * (vy[i, jdx_up] - vy[i, jdx_down]) * dy_scale_inv

        # (1/2)(huv)_x + (1/2)(-h_x*u*v + h*v_x*u - h*v*u_x)
        + half * (h[idx_right, j] * vx[idx_right, j] * vy[idx_right, j] - h[idx_left, j] * vx[idx_left, j] * vy[idx_left, j]) * dx_scale_inv  # (1/2)(huv)_x
        + half * (-(h[idx_right, j] - h[idx_left, j]) * dx_scale_inv * vx[i, j] * vy[i, j]                      # (1/2)(-h_x*u*v)
                  +
                  h[i, j] * (vy[idx_right, j] - vy[idx_left, j]) * dx_scale_inv * vx[i, j]                    # (1/2)(h*v_x*u) = (1/2)(h*u_x*v)
                  -
                  h[i, j] * vy[i, j] * (vx[idx_right, j] - vx[idx_left, j]) * dx_scale_inv)                   # (1/2)(-h*v*u_x)

        #dispersive terms
        + (λ / 6) * (eta[i, j]^2 / h[i, j]^2) * (h[i, jdx_up] - h[i, jdx_down]) * dy_scale_inv
        + (λ / 3) * (eta[i, jdx_up] - eta[i, jdx_down]) * dy_scale_inv
        -
        (λ / 3) * (eta[i, j] / h[i, j]) * (eta[i, jdx_up] - eta[i, jdx_down]) * dy_scale_inv
        -
        (λ / 6) * ((eta[i, jdx_up]^2 / h[i, jdx_up]) - (eta[i, jdx_down]^2 / h[i, jdx_down])) * dy_scale_inv

        # Bathymetrie (y-Richtung) 
        +
        (λ / 2) * (1 - eta[i, j] / h[i, j]) * b_y[i, j]
    )

    #  η-Gleichung 
    @inbounds deta[i, j] = -((eta[idx_right, j] - eta[idx_left, j]) * dx_scale_inv * vx[i, j] +
                             (eta[i, jdx_up] - eta[i, jdx_down]) * dy_scale_inv * vy[i, j] +
                             3 * half * b_x[i, j] * vx[i, j] +
                             3 * half * b_y[i, j] * vy[i, j]) + w[i, j]

    #  w-Gleichung 
    @inbounds dw[i, j] = minus_inv_h_ij * (
        # x-Terme
        half * ((h[idx_right, j] * vx[idx_right, j] * w[idx_right, j] - h[idx_left, j] * vx[idx_left, j] * w[idx_left, j]) * dx_scale_inv)
        -
        half * (h[idx_right, j] - h[idx_left, j]) * dx_scale_inv * vx[i, j] * w[i, j]
        -
        half * h[i, j] * (vx[idx_right, j] - vx[idx_left, j]) * dx_scale_inv * w[i, j]
        + half * h[i, j] * vx[i, j] * (w[idx_right, j] - w[idx_left, j]) * dx_scale_inv

        # y-Terme
        + half * ((h[i, jdx_up] * vy[i, jdx_up] * w[i, jdx_up] - h[i, jdx_down] * vy[i, jdx_down] * w[i, jdx_down]) * dy_scale_inv)
        -
        half * (h[i, jdx_up] - h[i, jdx_down]) * dy_scale_inv * vy[i, j] * w[i, j]
        -
        half * h[i, j] * (vy[i, jdx_up] - vy[i, jdx_down]) * dy_scale_inv * w[i, j]
        +
        half * h[i, j] * vy[i, j] * (w[i, jdx_up] - w[i, jdx_down]) * dy_scale_inv
    ) + λ / h[i, j] * (1 - eta[i, j] / h[i, j])

end


# reflecting boundary conditions
function get_grid(N, min, max, ::Val{true})
    return range(min, max, length=N)
end

# periodic boundary conditions
function get_grid(N, min, max, ::Val{false})
    full_range = range(min, max, length=N + 1)
    return range(full_range[1], stop=full_range[end-1],
        length=N)
end


# reflecting boundary conditions
function calculate_f_x_and_f_y(f, dx, dy, ::Val{true})

    Nx, Ny = size(f)
    f_x = zeros(Nx, Ny)
    f_y = zeros(Nx, Ny)
    half = eltype(f)(0.5)

    for j in 1:Ny
        jdx_down = (j == 1) ? 1 : j - 1 # use backward difference at bottom boundary
        jdx_up = (j == Ny) ? Ny : j + 1 # use forward difference at top boundary
        # Define scaling factors for finite differences
        dy_scale_inv = (j == 1 || j == Ny) ? 1 / dy : half / dy
        for i in 1:Nx
            idx_left = (i == 1) ? 1 : i - 1 # use backward difference at left boundary
            idx_right = (i == Nx) ? Nx : i + 1 # use forward difference at right boundary
            # Define scaling factors for finite differences
            dx_scale_inv = (i == 1 || i == Nx) ? 1 / dx : half / dx

            # Calculate derivatives
            f_x[i, j] = (f[idx_right, j] - f[idx_left, j]) * dx_scale_inv
            f_y[i, j] = (f[i, jdx_up] - f[i, jdx_down]) * dy_scale_inv
        end
    end

    return f_x, f_y
end

# periodic boundary conditions
function calculate_f_x_and_f_y(f, dx, dy, ::Val{false})

    Nx, Ny = size(f)
    f_x = zeros(Nx, Ny)
    f_y = zeros(Nx, Ny)
    half = eltype(f)(0.5)

    dx_scale_inv = half / dx
    dy_scale_inv = half / dy

    for j in 1:Ny
        jdx_down = mod1(j - 1, Ny)
        jdx_up = mod1(j + 1, Ny)
        for i in 1:Nx
            idx_left = mod1(i - 1, Nx)
            idx_right = mod1(i + 1, Nx)

            # Calculate derivatives
            f_x[i, j] = (f[idx_right, j] - f[idx_left, j]) * dx_scale_inv
            f_y[i, j] = (f[i, jdx_up] - f[i, jdx_down]) * dy_scale_inv
        end
    end

    return f_x, f_y
end

function u_analytic(t, coord, A, coord0, coordmin, coordmax, h∞=1.0)
    g = 9.81
    ϵ = A / h∞
    κ = sqrt(3 * ϵ / (4 * h∞^2 * (1 + ϵ)))
    c = sqrt(g * h∞ * (1 + ϵ))

    # Periodic boundary conditions
    L = coordmax - coordmin
    coord_t = coord - c * t
    coord_t = coordmin + mod(coord_t - coordmin, L)   # Apply periodicity
    coord_t -= coord0

    # coord_t = mod(coord - c * t - coordmin, coordmax - coordmin) + coordmin

    h = h∞ * (1 + ϵ * sech(κ * coord_t)^2)
    v = c * (1 - h∞ / h)

    return h, v
end

function setup_solitary_wave_2D(t, gridx, gridy, backend, reflecting_bc; coord0=0.0, A=0.3, h∞=1.0, direction=:x, b=zeros(length(gridx), length(gridy)))
    nx, ny = length(gridx), length(gridy)
    dx, dy = step(gridx), step(gridy)

    if direction == :x
        # Wave propagates in x direction
        coordmin, coordmax = extrema(gridx)
        if reflecting_bc == Val(false)
            # coordmax should be == xmax 
            coordmax = coordmax + step(gridx)
        end

        # Get 1D profile along x
        h_profile = zero(gridx)
        v_profile = zero(gridx)

        for (i, coord) in enumerate(gridx)
            h_profile[i], v_profile[i] = u_analytic(t, coord, A, coord0, coordmin, coordmax, h∞)
        end

        # Create 2D arrays (nx, ny) - constant in y direction
        h = repeat(h_profile, 1, ny)     # Each column is the same h profile
        vx = repeat(v_profile, 1, ny)    # Each column is the same vx profile  
        vy = zeros(nx, ny)               # No y-component velocity

    elseif direction == :y
        # Wave propagates in y direction
        coordmin, coordmax = extrema(gridy)
        if reflecting_bc == Val(false)
            # coordmax should be == ymax
            coordmax = coordmax + step(gridy)
        end

        # Get 1D profile along y
        h_profile = zero(gridy)
        v_profile = zero(gridy)

        for (i, coord) in enumerate(gridy)
            h_profile[i], v_profile[i] = u_analytic(t, coord, A, coord0, coordmin, coordmax, h∞)
        end

        # Create 2D arrays (nx, ny) - constant in x direction
        h = repeat(h_profile', nx, 1)    # Each row is the same h profile
        vx = zeros(nx, ny)               # No x-component velocity
        vy = repeat(v_profile', nx, 1)   # Each row is the same vy profile

    else
        error("Direction must be :x or :y")
    end

    h = h .- b  # adjust for bathymetry

    b_x, b_y = calculate_f_x_and_f_y(b, dx, dy, reflecting_bc)
    vx_x, _ = calculate_f_x_and_f_y(vx, dx, dy, reflecting_bc)
    _, vy_y = calculate_f_x_and_f_y(vy, dx, dy, reflecting_bc)

    w = @. -h * (vx_x + vy_y) + 1.5 * (vx * b_x + vy * b_y)

    q0 = zeros(nx, ny, 5)
    q0[:, :, 1] .= h
    q0[:, :, 2] .= vx
    q0[:, :, 3] .= vy
    q0[:, :, 4] .= h
    q0[:, :, 5] .= w

    return adapt(backend, q0)
end


function create_cache(; backend, λ, g, gridx, gridy, b, reflecting_bc, source_kernel=nothing)
    dx = step(gridx)
    dy = step(gridy)

    dx_inv = 1 / dx
    dy_inv = 1 / dy

    Nx, Ny = length(gridx), length(gridy)
    b_x, b_y = calculate_f_x_and_f_y(b, dx, dy, reflecting_bc)

    b = adapt(backend, b)
    b_x = adapt(backend, b_x)
    b_y = adapt(backend, b_y)
    gridx = adapt(backend, gridx)
    gridy = adapt(backend, gridy)

    return (; λ, g, Nx, Ny, dx, dx_inv, dy, dy_inv, b, b_x, b_y, gridx, gridy, reflecting_bc, source_kernel)
end


# callbacks
function save_func_all(u, t, integrator)
    return copy(adapt(CPU(), u))
end

function save_func_h(u, t, integrator)
    return copy(adapt(CPU(), u[:, :, 1]))
end

function save_and_print_callback(saveat; print_every_n=100, save_everything=false)
    # reset counter
    step_counter = Ref(0)
    # Callback that increments counter and prints every 100 steps
    function print_condition(u, t, integrator)
        step_counter[] += 1
        return step_counter[] % print_every_n == 0
    end

    function print_affect!(integrator)
        println("Step $(step_counter[]), t = $(integrator.t)")
    end

    print_cb = DiscreteCallback(print_condition, print_affect!, save_positions = (false,false))


    if save_everything
        saved_values = SavedValues(Float64, Array{Float64,3})
        save_cb = SavingCallback(save_func_all, saved_values, saveat=saveat)
    else
        saved_values = SavedValues(Float64, Matrix{Float64})
        save_cb = SavingCallback(save_func_h, saved_values, saveat=saveat)
    end



    return (CallbackSet(save_cb, print_cb), saved_values)

end

function gaussian(x, y; x0=0.0, y0=0.0, A=0.1, σ=1.0)
    return A * exp(-((x - x0)^2 + (y - y0)^2) / (2 * σ^2))
end




##########################################
# Function for analysis of the results
##########################################

function fit_powerlaw(N_values, errors)
    fit_function(N, params) = @. params[2] * N^(-params[1])
    fit = curve_fit(fit_function, N_values, errors, [2.0, 0.1])
    coef(fit)
    return coef(fit)
end

function calculate_2D_l2(u, v, gridx, gridy, reflecting_bc::Val{reflecting_bool}) where reflecting_bool
    dx, dy = step(gridx), step(gridy)

    x = @. (u - v)^2 * dx * dy

    if reflecting_bool
        x[begin, :] *= 0.5
        x[end, :] *= 0.5
        x[:, begin] *= 0.5
        x[:, end] *= 0.5
    end
    return sqrt(sum(x))
end

function calculate_1D_l2(u, v, gridx, gridy, directions)
    if directions == :x
        dx = step(gridx)
        # just look at one row/column
        res = @. (u[:, 1] - v[:, 1])^2 * dx
        return sqrt(sum(res))
    elseif directions == :y
        dy = step(gridy)
        # just look at one row/column
        res = @. (u[1, :] - v[1, :])^2 * dy
        return sqrt(sum(res))
    end
end

function calculate_2D_lmax(u, v)
    return maximum(abs, u .- v)
end

# numerically compute energy conservation

function compute_energy_partials_split_form(q, cache)
    h, vx, vy, eta, w = eachslice(q; dims=3)
    (; g, λ, b) = cache

    # ∂E/∂h
    dE_dh = @. 0.5 * (vx^2 + vy^2) +
               0.5 * g * (2 * h + 2 * b) +
               (1 / 6) * w^2 +
               (λ / 6) * (1 - eta^2 / h^2)

    # ∂E/∂vx
    dE_dvx = @. h * vx

    # ∂E/∂vy
    dE_dvy = @. h * vy

    # ∂E/∂eta
    dE_deta = @. -(λ / 3) * (1 - eta / h)

    # ∂E/∂w
    dE_dw = @. (1 / 3) * h * w

    return dE_dh, dE_dvx, dE_dvy, dE_deta, dE_dw
end





function compute_energy_time_derivative(q, cache, t)
    dE_dh, dE_dvx, dE_dvy, dE_deta, dE_dw = compute_energy_partials_split_form(q, cache)

    dq = similar(q)
    rhs_split!(dq, q, cache, t)  # calculate time derivatives

    dh, dvx, dvy, deta, dw = eachslice(dq; dims=3)

    (; dx, dy, reflecting_bc) = cache

    # dE/dt = Σ(∂E/∂qi * ∂qi/∂t)
    dE_dt_local = @. dE_dh * dh + dE_dvx * dvx + dE_dvy * dvy + dE_deta * deta + dE_dw * dw


    # integrate over the domain

    # for reflecting bc, weigh the boundaries with 1/2
    if reflecting_bc == Val(true)
        dE_dt_local[1, :] *= 0.5
        dE_dt_local[end, :] *= 0.5
        dE_dt_local[:, 1] *= 0.5
        dE_dt_local[:, end] *= 0.5
    end


    dE_dt = sum(dE_dt_local) * dx * dy

    return dE_dt
end