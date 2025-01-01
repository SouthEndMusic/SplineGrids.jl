"""
    RefinementMatrix(m, n, row_pointer, column_start, nzval)

The refinement matrix is a sparse matrix encoding for matrices which:

 1. Have consecutive non-zeros in all rows and columns
 2. Have at least 1 nonzero in every row and column

The nonzeros are stored in a dense vector per row.

## Fields

  - `m`: The number of rows of the matrix
  - `n`: The number of columns of the matrix
  - `row_pointer`: `row_pointer[i]` indicates where in `nzval` the data for the `i`-th row starts
  - `column_start`: `column_start[i]` indicates at which column index the first non-zero for the `i`-th row is
  - `nzval`: The nonzero values in the matrix
"""
struct RefinementMatrix{
    I <: AbstractVector{Ti} where {Ti},
    V <: AbstractVector{Tv} where {Tv},
    Tv,
    Ti <: Integer
} <: AbstractRefinementMatrix{Tv, Ti}
    m::Int
    n::Int
    row_pointer::I
    column_start::I
    nzval::V
    function RefinementMatrix(m, n, row_pointer, column_start, nzval)
        @assert length(row_pointer) == length(column_start) == m
        @assert row_pointer == sort(row_pointer)
        backend = get_backend(nzval)
        valid_row = KernelAbstractions.zeros(backend, Bool, length(row_pointer))
        backend = get_backend(nzval)
        validate_refinement_matrix_kernel(backend)(
            valid_row,
            row_pointer,
            column_start,
            length(nzval),
            ndrange = size(valid_row)
        )
        synchronize(backend)
        if !all(valid_row)
            error("Invalid rows: $(findall(x -> !x, valid_row)).")
        end
        new{
            typeof(row_pointer),
            typeof(nzval),
            eltype(nzval),
            eltype(row_pointer)
        }(
            m, n, row_pointer, column_start, nzval
        )
    end
end

# Get the index of the column with the last nonzero 
# for row i
function get_column_end(
        row_pointer::AbstractVector{<:Integer},
        column_start::Integer,
        n_nzval::Integer,
        i::Integer)
    row_pointerᵢ₊₁ = (i == length(row_pointer)) ? n_nzval + 1 :
                     row_pointer[i + 1]
    row_pointerᵢ = row_pointer[i]
    n_non_zeros = row_pointerᵢ₊₁ - row_pointerᵢ
    column_end = column_start + n_non_zeros - 1
    column_end
end

# Get the range of non-zeros start:end of column i
function get_column_range(
        row_pointer::AbstractVector{<:Integer},
        column_start::AbstractVector{<:Integer},
        n_nzval::Integer,
        i::Integer)
    column_start_ = column_start[i]
    column_end = get_column_end(row_pointer, column_start_, n_nzval, i)
    column_start_, column_end
end

# Validate:
# 1. That the non-zeros of the columns are consecutive
# 2. That every row and column has at least one nonzero
@kernel function validate_refinement_matrix_kernel(
        valid_row,
        @Const(row_pointer),
        @Const(column_start),
        n_nzval
)
    # Row index
    i = @index(Global, Linear)

    column_start_, column_end = get_column_range(
        row_pointer, column_start, n_nzval, i)

    # Validate that there is at least one nonzero in this row
    if column_end ≥ column_start_
        if (i == 1)
            valid_row[i] = true
        else
            column_start_prev, column_end_prev = get_column_range(
                row_pointer, column_start, n_nzval, i - 1
            )

            # Validate that the start of this row is between the start of the
            # previous row and the column after th end of the previous row
            condition_left_boundary = (column_start_prev ≤ column_start_ ≤
                                       column_end_prev + 1)

            # Validate that the end of this row is at least at the end of
            # the previous row
            if condition_left_boundary
                condition_right_boundary = (column_end ≥ column_end_prev)
                if condition_right_boundary
                    valid_row[i] = true
                end
            end
        end
    end
end

# Compute the non-zeros of the refinement matrix product C = A * B
@kernel function refinement_matrix_multiplication_kernel(
        nzval_C,
        @Const(row_pointer_C),
        @Const(column_start_C),
        @Const(row_pointer_A),
        @Const(row_pointer_B),
        @Const(column_start_A),
        @Const(column_start_B),
        @Const(nzval_A),
        @Const(nzval_B)
)
    # i: Row index (of A and C)
    i = @index(Global, Linear)

    # Source columns of A
    A_column_start, A_column_end = get_column_range(
        row_pointer_A, column_start_A, length(nzval_A), i)

    # Target columns of product C
    C_column_start, C_column_end = get_column_range(
        row_pointer_C, column_start_C, length(nzval_C), i)

    nzval_pointer_A = row_pointer_A[i]

    # k: Column index of A, row index of B
    for k in A_column_start:A_column_end
        B_column_start, B_column_end = get_column_range(
            row_pointer_B, column_start_B, length(nzval_B), k)

        nzval_pointer_C = row_pointer_C[i]

        # j: column index (of B and C)
        for j in C_column_start:C_column_end
            if B_column_start ≤ j ≤ B_column_end
                nzval_pointer_B = row_pointer_B[k] + j - B_column_start
                nzval_C[nzval_pointer_C] += nzval_A[nzval_pointer_A] *
                                            nzval_B[nzval_pointer_B]
            end
            # TODO: An optimization here is possible:
            # When a zero after a nonzero is found, the loop over j can be cut of 
            nzval_pointer_C += 1
        end
        nzval_pointer_A += 1
    end
end

# Compute for the refinement matrix product C = A * B
# per row the number of non-zeros and the column of the first nonzero
@kernel function refinement_matrix_mult_nonzeros_kernel(
        n_nonzero_C,
        column_start_C,
        @Const(row_pointer_A),
        @Const(row_pointer_B),
        @Const(column_start_A),
        @Const(column_start_B),
        n_nzval_A,
        n_nzval_B,
        n_columns_B
)
    # i: Row index (of A and C)
    i = @index(Global, Linear)

    A_column_start, A_column_end = get_column_range(
        row_pointer_A, column_start_A, n_nzval_A, i)

    n_non_zeros = 0
    column_start = 0

    # j: Column index (of B and C)
    for j in n_columns_B:-1:1
        # If there is overlap between the non-zeros of row i of matrix A 
        # and column j of matrix B then this yields a non-zero in row i of matrix C
        # k: Row index of B, column index of A
        for k in A_column_start:A_column_end
            B_column_start, B_column_end = get_column_range(
                row_pointer_B, column_start_B, n_nzval_B, k)

            if B_column_start ≤ j ≤ B_column_end
                n_non_zeros += 1
                column_start = j
                break
            end
        end
    end

    n_nonzero_C[i] = n_non_zeros
    column_start_C[i] = column_start
end

function Base.:*(
        A::AbstractRefinementMatrix{Tv, Ti},
        B::AbstractRefinementMatrix{Tv, Ti}
)::AbstractRefinementMatrix{Tv, Ti} where {Tv, Ti}
    if A.n != B.m
        throw(DimensionMismatch("Inner dimensions must match"))
    end

    backend = get_backend(A.nzval)

    n_nonzero_C = allocate(backend, Ti, A.m)
    column_start_C = allocate(backend, Ti, A.m)

    AB_args = (
        A.row_pointer,
        B.row_pointer,
        A.column_start,
        B.column_start
    )

    refinement_matrix_mult_nonzeros_kernel(backend)(
        n_nonzero_C,
        column_start_C,
        AB_args...,
        length(A.nzval),
        length(B.nzval),
        B.n,
        ndrange = (A.m,)
    )
    synchronize(backend)

    row_pointer_C = KernelAbstractions.zeros(backend, Ti, A.m)
    cumsum!(view(row_pointer_C, 2:(A.m)), n_nonzero_C[1:(end - 1)])
    row_pointer_C .+= 1
    nzval_C = KernelAbstractions.zeros(backend, Tv, sum(n_nonzero_C))

    refinement_matrix_multiplication_kernel(backend)(
        nzval_C,
        row_pointer_C,
        column_start_C,
        AB_args...,
        A.nzval,
        B.nzval,
        ndrange = (A.m,)
    )
    synchronize(backend)

    RefinementMatrix(
        A.m,
        B.n,
        row_pointer_C,
        column_start_C,
        nzval_C
    )
end

@kernel function collect_refinement_matrix_kernel(
        out,
        @Const(row_pointer),
        @Const(column_start),
        @Const(nzval)
)
    # Row index
    i = @index(Global, Linear)

    idx_data_start = row_pointer[i]
    idx_data_end = (i == length(row_pointer)) ? length(nzval) : row_pointer[i + 1] - 1

    idx_column = column_start[i]

    for idx_data in idx_data_start:idx_data_end
        out[i, idx_column] = nzval[idx_data]
        idx_column += 1
    end
end

function Base.collect(A::AbstractRefinementMatrix{Tv}) where {Tv}
    backend = get_backend(A.nzval)
    out = KernelAbstractions.zeros(backend, Tv, A.m, A.n)

    collect_refinement_matrix_kernel(backend)(
        out,
        A.row_pointer,
        A.column_start,
        A.nzval,
        ndrange = (A.m,)
    )
    synchronize(backend)

    out
end

function RefinementMatrix(A::Matrix{Tv}) where {Tv}
    row_pointer = Int[]
    column_start = Int[]
    nzval = Tv[]

    pointer = 1

    for row in eachrow(A)
        idx_first = findfirst(!iszero, row)
        idx_last = findlast(!iszero, row)
        push!(row_pointer, pointer)

        if isnothing(idx_first)
            push!(column_start, 0)
        else
            push!(column_start, idx_first)
            append!(nzval, view(row, idx_first:idx_last))
            pointer += idx_last - idx_first + 1
        end
    end

    RefinementMatrix(
        size(A)...,
        row_pointer,
        column_start,
        nzval
    )
end

Base.size(A::AbstractRefinementMatrix) = (A.m, A.n)
Base.length(A::AbstractRefinementMatrix) = A.m * A.n

function Base.:(==)(A::AbstractRefinementMatrix{Tv, Ti},
        B::AbstractRefinementMatrix{Tv, Ti}) where {Tv, Ti}
    (size(A) == size(B)) && (A.row_pointer == B.row_pointer) &&
        (A.column_start == B.column_start) && (A.nzval == B.nzval)
end

function Adapt.adapt(backend::Backend, A::AbstractRefinementMatrix)
    RefinementMatrix(
        size(A)...,
        adapt(backend, A.row_pointer),
        adapt(backend, A.column_start),
        adapt(backend, A.nzval)
    )
end