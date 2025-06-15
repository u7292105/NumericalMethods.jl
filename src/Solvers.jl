using ArgCheck

"""
    bisection(f, inita::Real, initb::Real; tol::Real = 1e-16, maxiter::Int = 1000)

Compute a root of the function `f` between `inita` and `initb` with tolerance
`tol` and maximum iterations `maxiter`. Might return `inita` or `initb` if 
tolerably close to a root.

Default values `tol = 1e-16` and `maxiter = 1000` are used if unspecified.
"""
function bisection(f, inita::Real, initb::Real; tol::Real = 1e-16, maxiter::Int = 1000)
    @argcheck inita < initb "lower bound must be smaller than upper bound"
    @argcheck tol > 0 "tolerance must be positive"
    a, b, c = inita, initb, nothing

    fa, fb, fc = f(a)::Real, f(b)::Real, nothing
    @argcheck fa * fb <= 0 "sign must change on the interval"
    if abs(fa) < tol return a end
    if abs(fb) < tol return b end

    for _ in 1:maxiter
        c = (a + b)/2
        fc = f(c)::Real

        if abs(fc) < tol
            return c
        elseif fa * fc < 0
            b, fb = c, fc
        else
            a, fa = c, fc
        end
    end
    error("bisection did not converge in $maxiter iterations, reached f($c) = $fc")
end