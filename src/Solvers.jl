using ArgCheck

"""
    function bisection(
        f, a::Real, b::Real;
        xtol::Real = 1e-16, ftol::Real = 1e-16, maxiter::Int = 1000
    )

Compute a root of the function `f` on the interval `[a, b]` with the bisection
algorithm.

# Arguments
- `f`: Real to Real function
- `a::Real`: lower interval bound
- `b::Real`: upper interval bound
- `xtol::Real = 1e-16`: domain tolerance
- `ftol::Real = 1e-16`: codomain tolerance
- `maxiter::Int = 1000`: maximum iterations

# Returns
- `(root::Real, f(root)::Real, noiter::Int, convergence::Bool)`
"""
function bisection(
    f, a::Real, b::Real;
    xtol::Real = 1e-16, ftol::Real = 1e-16, maxiter::Int = 1000
)
    @argcheck a < b "lower bound must be smaller than upper bound"
    @argcheck xtol > 0 && ftol > 0 "tolerances must be positive"
    @argcheck maxiter >= 1 "maximum iterations must be at least 1"
    u, v, e, n = f(a)::Real, f(b)::Real, b-a, 0
    @argcheck u * v <= 0 "sign must change on the interval"

    while true
        e = e/2
        m = a + e
        w = f(m)::Real
        n += 1

        if e < xtol
            return (m, w, n, true)
        elseif abs(w) < ftol
            return (m, w, n, true)
        elseif n == maxiter
            return (m, w, n, false)
        end

        if u * w <= 0
            b, v = m, w
        else
            a, u = m, w
        end
    end
end

"""
    function newton(
        f, df, x0::Real;
        xtol::Real = 1e-16, ftol::Real = 1e-16, maxiter = 1000
    )

Compute a root of the function `f` with initial guess `x0` with Newton's 
method.

# Arguments
- `f`: Real to Real function
- `df`: Real to Real derivative of `f`
- `x0::Real`: initial guess
- `xtol::Real = 1e-16`: domain tolerance
- `ftol::Real = 1e-16`: codomain tolerance
- `maxiter::Int = 1000`: maximum iterations

# Returns
- `(root::Real, f(root)::Real, noiter::Int, convergence::Bool)`
"""
function newton(
    f, df, x0::Real;
    xtol::Real = 1e-16, ftol::Real = 1e-16, maxiter::Int = 1000
)
    @argcheck xtol > 0 && ftol > 0 "tolerances must be positive"
    @argcheck maxiter >= 1 "maximum iterations must be at least 1"
    xold = x0
    fxold = f(x0)::Real
    dfxold = df(x0)::Real
    n = 0

    while true
        xnew = xold - (fxold / dfxold)
        fxnew = f(xnew)::Real
        dfxnew = df(xnew)::Real
        n += 1

        if abs(xnew - xold) < xtol
            return (xnew, fxnew, n, true)
        elseif abs(fxnew) < ftol
            return (xnew, fxnew, n, true)
        elseif n == maxiter
            return (xnew, fxnew, n, false)
        end
        xold, fxold, dfxold = xnew, fxnew, dfxnew
    end
end