@recipe function f(s::SplineDimension)
    data = decompress(s)
    (; knot_vector, sample_points) = s
    n_basis_functions = get_n_basis_functions(s)

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
        s.knot_vector.knot_values, zero(s.knot_vector.knot_values)
    end
end
