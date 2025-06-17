using Test
using NumericalMethods

@testset "bisection tests" begin
    poly(x) = x
    x, _, _, _, _ = bisection(poly, -1, 1)
    @test x == 0
end