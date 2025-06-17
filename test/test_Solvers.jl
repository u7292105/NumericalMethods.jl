using Test
using NumericalMethods

@testset "bisection" begin
    id(x) = x
    quad(x) = x*(x-2)
    frac(x) = (x^3 + 2x - 5) / (x^2 - 3)
    badf1(x) = x, x
    badf2(x, _) = x

    @testset "bisection errors" begin
        @test_throws ArgumentError bisection(id, 1, -1)
        @test_throws ArgumentError bisection(id, -1, 1; xtol=-1)
        @test_throws ArgumentError bisection(id, -1, 1; ftol=-1)
        @test_throws ArgumentError bisection(id, 1, 2)
        @test_throws TypeError bisection(badf1, -1, 1)
        @test_throws MethodError bisection(badf2, -1, 1)
    end

    @testset "bisection kwargs" begin
        xtolval = 1e-32
        ftolval = 1e-32
        x, _, _, _, _ = bisection(id, -1, 100; xtol=xtolval, ftol=ftolval)
        @test abs(x - 0) < 1e-32
        @test abs(id(x) - 0) < 1e-32

        maxiterval = 5
        x, _, _, n, conv_bool = bisection(id, -1, 100; maxiter=maxiterval)
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
            x, fx, _, n, conv_bool = bisection(f, a, b)
            @test abs(x - root) < 1e-16 || abs(f(x) - 0) < 1e-16
            @test abs(f(x) - fx) < 1e-16
            @test n > 0
            @test conv_bool
        end
    end
end