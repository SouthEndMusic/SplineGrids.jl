@kernel function expand_knot_vector_kernel(
        knots_all,
        @Const(knot_values),
        @Const(multiplicities)
)
    i = @index(Global, Linear)
    knot_value = knot_values[i]

    mult_sum = 0
    idx_end = i - 1
    for j in 1:idx_end
        mult_sum += multiplicities[j]
    end

    idx_start = mult_sum + 1
    idx_end = mult_sum + multiplicities[i]
    for k in idx_start:idx_end
        knots_all[k] = knot_value
    end
end

@kernel function set_sample_indices_kernel(
        sample_indices,
        @Const(sample_points),
        @Const(knots_all),
        degree
)
    i = @index(Global, Linear)
    sample_point = sample_points[i]

    sample_index = 0
    n_knots = length(knots_all)

    for knot in knots_all
        if sample_point < knot
            break
        else
            sample_index += 1
        end
    end

    sample_index = clamp(
        sample_index,
        degree + 1,
        n_knots - degree - 1
    )

    sample_indices[i] = sample_index
end

@kernel function decompress_basis_function_eval_kernel(
        decompressed_basis_functions,
        @Const(eval),
        @Const(sample_indices),
        degree,
        derivative_order
)
    i = @index(Global, Linear)
    l = sample_indices[i]

    idx_start = l - degree
    j = 0
    for k in idx_start:l
        j += 1
        decompressed_basis_functions[i, k] = eval[i, j, derivative_order + 1]
    end
end

@kernel function insert_kernel(out, @Const(v), i_insert, x)
    i = @index(Global, Linear)

    if i < i_insert
        out[i] = v[i]
    elseif i > i_insert
        out[i] = v[i - 1]
    else
        out[i] = x
    end
end