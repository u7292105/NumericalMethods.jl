using Test
using NumericalMethods

    id(x) = x
    quad(x) = x*(x-2)
    frac(x) = (x^3 + 2x - 5) / (x^2 - 3)

    did(x) = 1
    dquad(x) = 2x-2
    dfrac(x) = (x^4 - 11x^2 +10x - 6) / (x^2 - 3)^2

    badf1(x) = x, x
    badf2(x, _) = x

@testset "bisection" begin
    default_xtol = 1e-16
    default_ftol = 1e-16
    default_maxiter = 1000

    @testset "bisection errors" begin
        @test_throws ArgumentError bisection(id, 1, -1)
        @test_throws ArgumentError bisection(id, -1, 1; xtol=-1)
        @test_throws ArgumentError bisection(id, -1, 1; ftol=-1)
        @test_throws ArgumentError bisection(id, -1, 1; maxiter=-1)
        @test_throws ArgumentError bisection(id, 1, 2)
        @test_throws TypeError bisection(badf1, -1, 1)
        @test_throws MethodError bisection(badf2, -1, 1)
    end

    @testset "bisection kwargs" begin
        xtolval = 1e-32
        ftolval = 1e-32
        x, _, _, _ = bisection(id, -1, 100; xtol=xtolval, ftol=ftolval)
        @test abs(x - 0) < xtolval
        @test abs(id(x) - 0) < ftolval

        maxiterval = 5
        x, _, n, conv_bool = bisection(id, -1, 100; maxiter=maxiterval)
        @test n == maxiterval
        @test !conv_bool
    end

    @testset "bisection normal" begin
        test_args = [   (id, -1, 10, 0), 
                        (id, -10, 1, 0), 
                        (quad, -10, 1, 0), 
                        (frac, 1, 1.5, 1.32827)
                    ]
        for (f, a, b, root) in test_args
            x, fx, n, conv_bool = bisection(f, a, b)
            @test abs(x - root) < default_xtol || abs(f(x) - 0) < default_ftol
            @test abs(f(x) - fx) < 1e-16
            @test n > 0
            @test conv_bool
        end
    end
end

@testset "newton" begin
    default_xtol = 1e-16
    default_ftol = 1e-16
    default_maxiter = 1000

    @testset "newton errors" begin
        @test_throws ArgumentError newton(id, did, 1; xtol = -1)
        @test_throws ArgumentError newton(id, did, 1; ftol = -1)
        @test_throws ArgumentError newton(id, -1, 1; maxiter=-1)
        @test_throws TypeError newton(badf1, did, 1)
        @test_throws TypeError newton(id, badf1, 1)
        @test_throws MethodError newton(badf2, did, 1)
        @test_throws MethodError newton(id, badf2, 1)
    end

    @testset "newton kwargs" begin
        xtolval = 1e-32
        ftolval = 1e-32
        x, _, _, _ = newton(id, did, 1; xtol=xtolval, ftol=ftolval)
        @test abs(x - 0) < xtolval
        @test abs(id(x) - 0) < ftolval

        maxiterval = 5
        x, _, n, conv_bool = newton(quad, dquad, 1000; maxiter=maxiterval)
        @test n == maxiterval
        @test !conv_bool
    end

    @testset "newton normal" begin
        test_args = [   (id, did, 1, 0),
                        (id, did, 100, 0),
                        (id, did, -1, 0),
                        (id, did, -100, 0),
                        (quad, dquad, 0.5, 0),
                        (quad, dquad, 1.5, 2),
                        (quad, dquad, 100, 2),
                        (quad, dquad, -100, 0),
                        (frac, dfrac, 0, 1.32827),
                        (frac, dfrac, 2, 1.32827),
                        (frac, dfrac, -5, 1.32827),
                        (frac, dfrac, 5, 1.32827)
                    ]
        for (f, df, x0, root) in test_args
            x, fx, n, conv_bool = newton(f, df, x0)
            @test abs(x - root) < default_xtol || abs(f(x) - fx) < default_ftol
            @test abs(f(x) - fx) < 1e-16
            @test n > 0
            @test conv_bool
        end
    end
end