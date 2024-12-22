using SplineGrids
using Plots

n_control_points = (5, 6)
degree = (3, 2)
n_sample_points = (32, 55)

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)

dimensionalities = [
    ((1, Nout) for Nout in 1:4)...,
    ((2, Nout) for Nout in 1:3)...
]

for (Nin, Nout) in dimensionalities
    @testset "Nin = $Nin, Nout = $Nout" begin
        spline_grid = SplineGrid(spline_dimensions[1:Nin], Nout)
        @test try
            plot(spline_grid)
            true
        catch
            false
        end
    end
end