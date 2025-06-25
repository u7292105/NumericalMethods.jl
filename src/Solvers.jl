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

"""
    function secant(
        f, x0::Real;
        x1::Union{Nothing,Real} = nothing, xtol::Real = 1e-16, 
        ftol::Real = 1e-16, maxiter::Int = 1000, delta::Real = 1e-4
    )

Compute a root of the function `f` with initial guesses `x0` and `x1` with the 
Secant method. Second initial guess `x1` is optional.

# Arguments
- `f`: Real to Real function
- `x0::Real`: 1st initial guess
- `x1::Union{Nothing,Real} = nothing`: Optional 2nd initial guess
- `xtol::Real = 1e-16`: domain tolerance
- `ftol::Real = 1e-16`: codomain tolerance
- `maxiter::Int = 1000`: maximum iterations
- `delta::Real = 1e-4`: delta for estimating the 2nd initial guess
    - `x1` is set to `estimate_x1(f, x0, delta)` if `x1` is `nothing`

# Returns
- `(root::Real, f(root)::Real, noiter::Int, convergence::Bool)`
"""
function secant(
    f, x0::Real;
    x1::Union{Nothing,Real} = nothing, xtol::Real = 1e-16, 
    ftol::Real = 1e-16, maxiter::Int = 1000, delta::Real = 1e-4
)
    @argcheck xtol > 0 && ftol > 0 "tolerances must be positive"
    @argcheck maxiter >= 1 "maximum iterations must be at least 1"
    @argcheck delta > 0 "delta must be positive"
    if isnothing(x1)
        x1 = _estimate_x1(f, x0, delta)
    end
    xk0, xk1 = x0, x1
    fxk0, fxk1 = f(x0)::Real, f(x1)::Real
    n = 0

    while true
        xk2 = xk1 - fxk1 * (xk1 - xk0) / (fxk1 - fxk0)
        fxk2 = f(xk2)::Real
        n += 1

        if abs(xk2 - xk1) < xtol
            return (xk2, fxk2, n, true)
        elseif abs(fxk2 - fxk1) < ftol
            return (xk2, fxk2, n, true)
        elseif n == maxiter
            return (xk2, fxk2, n, false)
        end
        xk0, xk1 = xk1, xk2
        fxk0, fxk1 = fxk1, fxk2
    end
end

"""
    function _estimate_x1(f, x0::Real, delta::Real)

Estimate the 2nd initial guess for the Secant method. `x1` takes form `x0 + delta * abs(x0)` or `x0 - delta * abs(x0)` dependant on sign changes over the
intervals. `x1` defaults to the former if the sign changes on neither interval or both intervals.

# Arguments
- `f`: Real to Real function
- `x0::Real`: 1st initial guess
- `delta::Real`: delta for estimating the 2nd initial guess

# Returns
- `x1::Real`
"""
function _estimate_x1(f, x0::Real, delta::Real)
    x1pos = x0 + delta * abs(x0)
    x1neg = x0 - delta * abs(x0)
    if f(x0) * f(x1pos) < 0
        return x1pos
    elseif f(x0) * f(x1neg) < 0
        return x1neg
    else
        return x1pos
    end
end