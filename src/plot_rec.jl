const CONTROL_POINT_COLOR = :green
const CONTROL_POLYGON_COLOR = :orange
const SPLINE_GRID_COLOR = :blue

@recipe function f(spline_dimension::SplineDimension; derivative_order = 0)
    data = decompress(spline_dimension; derivative_order)
    (; sample_points, knot_vector) = spline_dimension
    (; knot_values) = knot_vector
    n_basis_functions = get_n_basis_functions(spline_dimension)

    # Plot each basis function
    for j in 1:n_basis_functions
        @series begin
            seriestype := :line
            label := "Basis function $j"
            sample_points, data[:, j]
        end
    end

    # Plot the sum of basis functions
    @series begin
        seriestype := :line
        label := "Basis function sum"
        sample_points, sum(data, dims = 2)[:]
    end

    # Plot the knots
    @series begin
        seriestype := :scatter
        label := "knots"
        knot_values, zero(knot_values)
    end
end

@recipe function f(spline_grid::AbstractSplineGrid{1, 1})
    (; eval, spline_dimensions) = spline_grid
    (; sample_points, knot_vector) = only(spline_dimensions)
    (; knot_values) = knot_vector

    xlabel --> "Input dimension 1"
    ylabel --> "Output dimension 1"

    # Plot the spline curve
    @series begin
        seriestype := :line
        label := "Spline curve"
        c := SPLINE_GRID_COLOR
        sample_points, eval
    end

    Min = minimum(eval)
    Max = maximum(eval)

    for (i, knot_value) in enumerate(knot_values)
        @series begin
            seriestype := :line
            label := (i == 1) ? "Knots" : ""
            color := :gray
            ls := :dash
            [knot_value, knot_value], [Min, Max]
        end
    end
end

@recipe function f(spline_grid::AbstractSplineGrid{1, 2})
    (; eval, control_points) = spline_grid

    xlabel --> "Output dimension 1"
    ylabel --> "Output dimension 2"

    # Plot the spline curve
    @series begin
        seriestype := :path
        label := "Spline curve"
        c := SPLINE_GRID_COLOR
        eval[:, 1], eval[:, 2]
    end

    # Plot the control polygon
    @series begin
        seriestype := :path
        label := "Control polygon"
        c := CONTROL_POLYGON_COLOR
        control_points[:, 1], control_points[:, 2]
    end

    # Plot the control points
    @series begin
        seriestype := :scatter
        label := "Control points"
        c := CONTROL_POINT_COLOR
        spline_grid.control_points[:, 1], spline_grid.control_points[:, 2]
    end
end

@recipe function f(spline_grid::AbstractSplineGrid{1, 3})
    (; eval, control_points) = spline_grid

    xlabel --> "Output dimension 1"
    ylabel --> "Output dimension 2"
    zlabel --> "Output dimension 3"

    # Plot the control polygon
    @series begin
        seriestype := :path
        label := "Control polygon"
        c := CONTROL_POLYGON_COLOR
        control_points[:, 1], control_points[:, 2], control_points[:, 3]
    end

    # Plot the control points
    @series begin
        seriestype := :scatter
        label := "Control points"
        c := CONTROL_POINT_COLOR
        control_points[:, 1], control_points[:, 2], control_points[:, 3]
    end

    # Plot the spline curve
    @series begin
        seriestype := :path
        label := "Spline curve"
        c := SPLINE_GRID_COLOR
        eval[:, 1], eval[:, 2], eval[:, 3]
    end

    return nothing
end

@recipe function f(spline_grid::AbstractSplineGrid{1, 4})
    (; eval, control_points) = spline_grid

    xlabel --> "Output dimension 1"
    ylabel --> "Output dimension 2"
    zlabel --> "Output dimension 3"

    # Plot the control polygon
    @series begin
        seriestype := :path
        label := "Control polygon"
        c := CONTROL_POLYGON_COLOR
        control_points[:, 1], control_points[:, 2], control_points[:, 3]
    end

    # Plot the control points
    @series begin
        seriestype := :scatter
        label := "Control points"
        c := CONTROL_POINT_COLOR
        control_points[:, 1], control_points[:, 2], control_points[:, 3]
    end

    # Plot the spline curve
    @series begin
        seriestype := :path
        label := "Spline curve"
        line_z := eval[:, 4]
        c := :viridis
        colorbar_title := "Output dimension 4"
        eval[:, 1], eval[:, 2], eval[:, 3]
    end

    return nothing
end

@recipe function f(spline_grid::AbstractSplineGrid{2, 1})
    (; spline_dimensions, eval) = spline_grid

    xlabel --> "Input dimension 1"
    ylabel --> "Input dimension 2"

    knot_values_1 = spline_dimensions[1].knot_vector.knot_values
    knot_values_2 = spline_dimensions[2].knot_vector.knot_values
    extent_1 = [first(knot_values_1), last(knot_values_1)]
    extent_2 = [first(knot_values_2), last(knot_values_2)]

    # Plot spline surface
    @series begin
        seriestype := :heatmap
        label := "Spline surface"
        colorbar_title := "Output dimension 1"
        spline_dimensions[1].sample_points, spline_dimensions[2].sample_points, eval[:, :]'
    end

    # Plot knot vectors
    for (i, knot_value) in enumerate(knot_values_1)
        @series begin
            seriestype := :line
            label := (i == 1) ? "Knots dimension 1" : ""
            color := :gray
            ls := :dash
            [knot_value, knot_value], extent_2
        end
    end

    for (i, knot_value) in enumerate(knot_values_2)
        @series begin
            seriestype := :line
            label := (i == 1) ? "Knots dimension 2" : ""
            color := :black
            ls := :dash
            extent_1, [knot_value, knot_value]
        end
    end
end

@recipe function f(spline_grid::AbstractSplineGrid{2, 2})
    (; eval, spline_dimensions) = spline_grid

    xlabel --> "Input dimension 1"
    ylabel --> "Input dimension 2"

    sample_points_1 = spline_dimensions[1].sample_points
    sample_points_2 = spline_dimensions[2].sample_points

    eval_u = @view eval[:, :, 1]
    eval_v = @view eval[:, :, 2]

    vec_length = min(minimum(diff(sample_points_1)), minimum(diff(sample_points_2)))

    eval_norms = @. sqrt(eval_u^2 + eval_v^2)

    eval_normalized = copy(eval)
    eval_normalized_u = @view eval_normalized[:, :, 1]
    eval_normalized_v = @view eval_normalized[:, :, 2]
    @. eval_normalized_u *= vec_length / eval_norms
    @. eval_normalized_v *= vec_length / eval_norms

    @series begin
        seriestype := :quiver
        label := "Spline vector field"
        quiver := (eval_normalized_u[:], eval_normalized_v[:])
        line_z := repeat(eval_norms[:], inner = 4)
        c := :viridis
        colorbar_title := "Vector norm"
        repeat(sample_points_1, 1, length(sample_points_2))[:],
        repeat(sample_points_2, 1, length(sample_points_1))'[:]
    end
end

@recipe function f(spline_grid::AbstractSplineGrid{2, 3})
    (; eval) = spline_grid

    xlabel --> "Output dimension 1"
    ylabel --> "Output dimension 2"
    zlabel --> "Output dimension 3"

    # Plot the spline surface
    @series begin
        seriestype := :surface
        colorbar := false
        eval[:, :, 1], eval[:, :, 2], eval[:, :, 3]
    end

    # TODO: For some reason combining anything with the surface plot yields an error in Plots.jl
    # # Plot the control points
    # @series begin
    #     seriestype := :scatter
    #     color := CONTROL_POINT_COLOR
    #     control_points[:, :, 1], control_points[:, :, 2], control_points[:, :, 3]
    # end

    # Plot the control polygon

    return nothing
end