# PDE solving example

In this section we solve the following PDE problem with a spline grid:

```math
\begin{align*}
\begin{cases}
    u : \Omega \rightarrow \mathbb{R} \\
    \Delta u = u^3 &\text{ for }& \mathbf{x} \in \Omega \\
    u = g(\mathbf{x}) \; &\text{ for }& \mathbf{x} \in \partial\Omega
\end{cases}
\end{align*},
```

where the domain is given by the open unit cube: $\Omega = (0,1)^3$. We assume that $g(x,y,z) = 0$ for $z \in (0,1)$ and $g(x,y,z) = f(x,y)$ for $z = 0, 1$.

We solve this problem by sampling the domain and enforcing the PDE on the interior points and the boundary condition on the boundary points. 

## The residual kernel

We first define the kernel which calculates the residual of a given approximation to the solution of the problem above. Note that this kernel is agnostic of the fact that the solution will come from a spline grid.

```@setup tutorial
using Base.Threads

# Set the number of threads to the number of available CPU cores
ENV["JULIA_NUM_THREADS"] = string(nthreads())
```

```@example tutorial
using KernelAbstractions

@kernel function pde_residual_kernel(
    residual,
    @Const(f),
    @Const(u),
    @Const(∂₁²u),
    @Const(∂₂²u),
    @Const(∂₃²u),
)
    I = @index(Global, Cartesian)

    is_boundary = false
    for (i, i_max) in zip(Tuple(I)[1:end-1], size(residual)[1:end-1])
        if (i == 1) || (i == i_max)
            is_boundary = true
            break
        end
    end

    residual[I] = if is_boundary
        if I[3] == 1 || I[3] == size(residual)[3]
            u[I] - f[I[1], I[2]]
        else
            u[I]
        end
    else
        ∂₁²u[I] + ∂₂²u[I] + ∂₃²u[I] - u[I]^3
    end
end
```

## The residual as a function of the spline grid control points

Whe now define a function which generates the input of the above residual kernel from a spline grid and computes the residual in place. It assumes a parameter object with the spline grid, `f` (an array that specifies the boundary condition for the bottom and top of the domain, the boundary condition is `0` elsewhere) and caches for the spline grid evaluation.

```@example tutorial
function pde_residual!(residual, control_points_flat, p)::Nothing
    (; spline_grid, f, u, ∂₁²u, ∂₂²u, ∂₃²u) = p

    control_points = reshape(
        control_points_flat,
        size(spline_grid.control_points)
    )

    evaluate!(spline_grid; control_points, eval=u)
    evaluate!(spline_grid; control_points, eval=∂₁²u, derivative_order=(2, 0, 0))
    evaluate!(spline_grid; control_points, eval=∂₂²u, derivative_order=(0, 2, 0))
    evaluate!(spline_grid; control_points, eval=∂₃²u, derivative_order=(0, 0, 2))

    pde_residual_kernel(get_backend(u))(
        reshape(residual, size(u)),
        f,
        u,
        ∂₁²u,
        ∂₂²u,
        ∂₃²u,
        ndrange=size(u)
    )

    return nothing
end
```

## Defining the spline grid

We would like to define a residual function $\mathbb{R}^N \rightarrow \mathbb{R}^N$, i.e. the number of sample points must equal the number of control points. For each dimension the maximum derivative order is 2, to evaluate the result of the Laplacian operator on the spline.

```@example tutorial
using SplineGrids

n_control_points = (25, 25, 25)
degree = (2, 2, 2)
n_sample_points = n_control_points
dim_out = 1

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points; max_derivative_order = 2)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid
```

## Defining the boundary condition

We define a boundary condition inspired by the Julia logo.

```@example tutorial
using GLMakie

function hill(x, y, x0, y0, β)
    r = sqrt((x - x0)^2 + (y - y0)^2)
    return exp(-(r / β)^2)
end

function hills(x, y)
    R = 0.25
    out = 0.0
    for θ in range(0, 2π, length=4)[1:end-1]
        x0 = 0.5 + R * cos(θ)
        y0 = 0.5 + R * sin(θ)
        out += hill(x, y, x0, y0, R / 2)
    end
    return 100 * out
end


f = [
    hills.(x, y) for
    x = spline_dimensions[1].sample_points,
    y = spline_dimensions[2].sample_points
]

heatmap(f)
```

## Defining the problem

We will solve the problem with a Jacobian free Newton-Krylov method. To do this, we need
to provide a Jacobian-vector product (JVP) function.

```@example tutorial
using NonlinearSolve
using Enzyme

p = (;
    spline_grid,
    f,
    u = similar(spline_grid.eval),
    ∂₁²u = similar(spline_grid.eval),
    ∂₂²u = similar(spline_grid.eval),
    ∂₃²u = similar(spline_grid.eval)
)

meta_p = (; 
    residual = similar(spline_grid.eval),
    p_duplicated = DuplicatedNoNeed(p, make_zero(p))
)

function pde_residual_jvp!(Jv, v, control_points_flat, meta_p)::Nothing
    for val in values(meta_p.p_duplicated.dval)
        make_zero!(val)
    end
    autodiff(
        Forward,
        pde_residual!,
        DuplicatedNoNeed(vec(meta_p.residual), Jv),
        Duplicated(control_points_flat, v),
        meta_p.p_duplicated
    )
    return nothing
end

nonlinear_function = NonlinearFunction(
    (residual, control_points_flat, p_) -> pde_residual!(residual, control_points_flat, p);
    jvp = pde_residual_jvp!
)

# Initial guess
x0 = zeros(length(spline_grid.control_points))

problem = NonlinearProblem(
    nonlinear_function,
    x0,
    meta_p
)
```

## Solving the problem

```@example tutorial
sol = solve(problem, NewtonRaphson(linsolve = KrylovJL_GMRES()))
```

## Viewing the solution

```@example tutorial
spline_grid.control_points .= reshape(sol.u, size(spline_grid.control_points))
evaluate!(spline_grid)
fig = Figure()
ax, plt = volume(fig[1,1], log.(spline_grid.eval[:, :, :, 1] .+ 1))
Colorbar(fig[1,2], plt; label = L"\log(u + 1)")
fig
```