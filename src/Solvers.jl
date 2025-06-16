using ArgCheck

"""
    function bisection(
        f, a::Real, b::Real;
        xtol::Real = 1e-16, ftol::Real = 1e-16, maxiter::Int = 1000
    )

Compute a root of the function `f` on the interval `[a, b]`.

# Arguments
- `f`: Real to Real function
- `a::Real`: lower interval bound
- `b::Real`: upper interval bound
- `xtol::Real = 1e-16`: domain tolerance
- `ftol::Real = 1e-16`: codomain tolerance
- `maxiter::Int = 1000`: maximum iterations

# Returns
- `(root::Real, f(root)::Real, xerror::Real, noiter::Int, convergence::Bool)`
"""
function bisection(
    f, a::Real, b::Real;
    xtol::Real = 1e-16, ftol::Real = 1e-16, maxiter::Int = 1000
)
    @argcheck a < b "lower bound must be smaller than upper bound"
    @argcheck xtol > 0 && ftol > 0 "tolerances must be positive"
    u, v, e, n = f(a)::Real, f(b)::Real, b-a, 0
    @argcheck u * v <= 0 "sign must change on the interval"

    while true
        e = e/2
        m = a + e
        w = f(m)::Real
        n += 1

        if e < xtol
            return (m, w, e, n, true)
        elseif abs(w) < ftol
            return (m, w, e, n, true)
        elseif n == maxiter
            return (m, w, e, n, false)
        end

        if u * w <= 0
            b, v = m, w
        else
            a, u = m, w
        end
    end
end