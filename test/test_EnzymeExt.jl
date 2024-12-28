using SplineGrids
using Enzyme
using Random

@testset "Enzyme" begin
    Random.seed!(1)

    n_control_points = (10, 10)
    degree = (2, 2)
    n_sample_points = (50, 50)
    Nout = 2

    spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
    spline_grid = SplineGrid(spline_dimensions, Nout)

    function loss(control_points_flat, spline_grid)
        evaluate!(spline_grid;
            control_points = reshape(control_points_flat, size(spline_grid.control_points))
        )
        return sum(spline_grid.eval)
    end

    control_points_flat = rand(length(spline_grid.control_points))

    dcontrol_points_flat = Duplicated(control_points_flat, make_zero(control_points_flat))
    dspline_grid = Duplicated(spline_grid, make_zero(spline_grid))

    @test try
        autodiff(
            Reverse,
            loss,
            Active,
            dcontrol_points_flat,
            dspline_grid
        )
        true
    catch
        false
    end

    # Test that a nonzero gradient was computed
    @test any(!iszero, dcontrol_points_flat.dval)

    @test try
        make_zero!(dspline_grid.dval)
        true
    catch
        false
    end
end