using SplineGrids
using KernelAbstractions
using Adapt

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

@testset "Nin = 1, Nout = 1" begin
    n_sample_points = 100
    for n_basis_functions in 2:10
        for degree in 1:(n_basis_functions - 1)
            spline_dimension = SplineDimension(
                n_basis_functions,
                degree,
                n_sample_points;
                backend
            )
            spline_grid = SplineGrid(spline_dimension, 1)
            spline_grid.control_points .= 1
            evaluate!(spline_grid)
            @test all(spline_grid.eval .≈ 1)
        end
    end
end

@testset "Nin = 2, Nout = 2" begin
    n_control_points = (5, 6)
    degree = (3, 2)
    n_sample_points = (7, 9)
    Nout = 2

    spline_dimensions = SplineDimension.(
        n_control_points,
        degree,
        n_sample_points;
        backend
    )
    spline_grid = SplineGrid(spline_dimensions, Nout)
    copyto!(
        spline_grid.control_points,
        [0.1286512149446628 0.989330939952115 0.733492448545175 0.5794756343618787 0.11070882026487838 0.26133752996006376; 0.11113515457625123 0.23135826137623983 0.46967245778179434 0.1236183748627232 0.8169302110385438 0.8941133259548234; 0.3310189719463331 0.5728511097513352 0.543496916829333 0.49787004423235026 0.27541446537092906 0.044945104057353524; 0.3663497772286072 0.010888588395355225 0.9031007598650994 0.38269754883935825 0.3301652647427398 0.22538336730505104; 0.11121593400768026 0.8498281072602938 0.15141863475808437 0.981856057718664 0.03554914225818451 0.6044665027464156;;;
         0.6372026128560905 0.22093519241233173 0.05608637192615418 0.8968022994859238 0.5826027419132229 0.6943761278711217; 0.8873259697643807 0.34896625070444054 0.6677168115009738 0.019325455816958326 0.8442963293878447 0.09315653801718138; 0.0751997400178176 0.03404959543348307 0.5865499242390589 0.5849387336372873 0.7504385341534244 0.7030718228624017; 0.6275159834003844 0.6265379095984482 0.650357679885941 0.06456851457005197 0.24252029490114768 0.6443206523975683; 0.5628895658150487 0.5805570340970917 0.9761588549617932 0.23835489276149868 0.16225413295209445 0.11017417901314286]
    )
    evaluate!(spline_grid)
    @test spline_grid.eval ≈
          adapt(backend,
        [0.1286512149446628 0.7421811972743845 0.8614116942486449 0.7462201581981305 0.6564840414535269 0.5401318843726657 0.3450922273133785 0.20696184945079976 0.26133752996006376; 0.14719165455423897 0.41833156712267155 0.5297966852821444 0.5220394475964073 0.4355122926292099 0.3690001675707926 0.421288019611396 0.5193986587676646 0.5903548950562434; 0.21212713529380467 0.3513801723101711 0.4542354151673586 0.4876542989307811 0.41859825866585276 0.3603595140655854 0.4262302848229905 0.5102941988400972 0.5066348840189345; 0.28488071892438116 0.36495494215849794 0.4809645150724781 0.5515189809148251 0.4952278829340427 0.41156336033841856 0.39999755233624046 0.38782662023536413 0.3023467253436454; 0.32491010532681536 0.3364896706929626 0.4881970244775462 0.6193236807358246 0.5691611535230561 0.45189057815870876 0.38169308996225093 0.31847185754969376 0.2221300495370486; 0.28381154686125665 0.3586416107536558 0.4820516161100149 0.5857022594759121 0.6012542373969256 0.5383217804480385 0.40651911920423295 0.3018969814951609 0.320506095150474; 0.11121593400768026 0.5778738798843641 0.5006233710091891 0.342524496690933 0.5666373462383743 0.7597630154160316 0.5087025999884243 0.2960668468128022 0.6044665027464156;;;
         0.6372026128560905 0.3043959449624992 0.13851078216924295 0.18178196543189756 0.476444335706039 0.7524378638443651 0.7397025206995733 0.6498210330992852 0.6943761278711217; 0.7055340750988145 0.40634653324413955 0.37429465610306395 0.435087515051254 0.41443418146437594 0.41789698225997274 0.5510382443555875 0.6011344645547687 0.35546213966106627; 0.5580294650569024 0.3713152769778062 0.4309327927946016 0.5243885904298894 0.43918924780627 0.3696194995296285 0.5099640802058502 0.6006602956337773 0.38214545161225244; 0.4163103583001001 0.34498981131609713 0.441847211379361 0.5388881508756761 0.46811822219082716 0.3937967705759495 0.4801831412821783 0.5774837991537468 0.5359052090348883; 0.43980602493064325 0.4304824338245724 0.5124558419069737 0.5607104814884867 0.4502305848797504 0.33432745170893385 0.36631238160420587 0.47906630966951136 0.6054701710087955; 0.5391765131807299 0.5606040917206397 0.636160258083579 0.6436651705361933 0.4609389873451282 0.2659394047881031 0.23662411914283762 0.32637340868549797 0.4885675516922504; 0.5628895658150487 0.6255903946346686 0.7783579445294424 0.8344831320786688 0.607256873861646 0.32106779306036 0.20030451285679657 0.1587467394435321 0.11017417901314286]
    )
end

@testset "Derivative order validation" begin
    using Logging

    n_control_points = (5, 6)
    degree = (3, 2)
    n_sample_points = (7, 9)
    Nout = 2
    max_derivative_order = 1

    spline_dimensions = SplineDimension.(
        n_control_points,
        degree,
        n_sample_points;
        max_derivative_order,
        backend
    )
    spline_grid = SplineGrid(spline_dimensions, Nout)

    logger = TestLogger()
    with_logger(logger) do
        @test_throws "Invalid derivative order(s) supplied. If you want to evaluate (higher order) derivatives, specify this at construction as SplineDimension(...; max_derivative_order)." evaluate!(
            spline_grid; derivative_order = (2, 1))
    end

    @test length(logger.logs) == 1
    @test logger.logs[1].level == Error
    @test logger.logs[1].message ==
          "The maximum derivative order available for spline dimension 1 is 1, got 2."
end