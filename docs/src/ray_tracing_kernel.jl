using StaticArraysCore
using LinearAlgebra

# Kernel function for computing the contribution of a ray intersection to the
# pixels close to the intersection.
function F(x, x0, w)
    if x < x0 - w
        -one(x) / 2
    elseif x > x0 + w
        one(x) / 2
    else
        x_transformed = (x - x0) / w
        (sin(π * x_transformed) / π + x_transformed) / 2
    end
end

@kernel function ray_tracing_kernel(
        render,
        @Const(u),
        @Const(∂₁u),
        @Const(∂₂u),
        @Const(x),
        @Const(y),
        scene_params
)
    (; r,
    z_screen,
    screen_size,
    ray_kernel_size) = scene_params

    T = eltype(u)
    TVec = SVector{3, T}
    I = @index(Global, Cartesian)

    # Origin of the ray on the lens surface
    o = TVec(x[I[1]], y[I[2]], u[I])

    # Lens surface normal at intersection
    n = normalize(TVec(∂₁u[I], ∂₂u[I], -1))

    # incoming light ray direction
    l = TVec(zero(T), zero(T), one(T))

    c = -n ⋅ l
    sqrt_arg = 1 - r^2 * (1 - c^2)

    if sqrt_arg >= 0
        # Refracted ray direction
        v = r * l + (r * c - √sqrt_arg) * n

        # Screen intersection 'time' for unit speed
        t_screen_int = (z_screen - o[3]) / v[3]

        if t_screen_int >= 0
            # Screen intersection point
            s = o + t_screen_int * v

            # Pixel size
            w_screen, h_screen = screen_size
            n_x, n_y = size(render)
            w_pixel = w_screen / n_x
            h_pixel = h_screen / n_y

            # Pixel intersection indices
            n_x, n_y = size(render)

            i = 1 + Int(floor((w_screen / 2 + s[1]) / w_pixel))
            j = 1 + Int(floor((h_screen / 2 + s[2]) / h_pixel))

            # Render contribution from this ray
            i_min = max(i - ray_kernel_size[1] - 1, 1)
            i_max = min(i + ray_kernel_size[1] + 1, n_x)
            j_min = max(j - ray_kernel_size[2] - 1, 1)
            j_max = min(j + ray_kernel_size[2] + 1, n_y)

            w_kernel = (ray_kernel_size[1] + 0.5) * w_pixel
            h_kernel = (ray_kernel_size[2] + 0.5) * h_pixel

            for i_ in i_min:i_max
                contribution_x = F(-0.5w_screen + i_ * w_pixel, s[1], w_kernel) -
                                 F(-0.5w_screen + (i_ - 1) * w_pixel, s[1], w_kernel)

                for j_ in j_min:j_max
                    contribution_y = F(-0.5h_screen + j_ * h_pixel, s[2], h_kernel) -
                                     F(
                        -0.5h_screen + (j_ - 1) * h_pixel, s[2], h_kernel)

                    Atomix.@atomic render[i_, j_] += contribution_x * contribution_y
                end
            end
        end
    end
end
