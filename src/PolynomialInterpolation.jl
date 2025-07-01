using ArgCheck

struct BasisPolynomial{F}
    func::F
    func_str::String
end
(bp::BasisPolynomial)(x::Real) = bp.func(x)

abstract type InterpolationPolynomial end

struct StandardPolynomial <: InterpolationPolynomial
    coeffs::Vector{Real}
    base_polys::Vector{BasisPolynomial}
end
struct LagrangePolynomial <: InterpolationPolynomial
    coeffs::Vector{Real}
    base_polys::Vector{BasisPolynomial}
end
struct NewtonPolynomial <: InterpolationPolynomial
    coeffs::Vector{Real}
    base_polys::Vector{BasisPolynomial}
end

function (ip::InterpolationPolynomial)(x::Real)
    return sum(ip.coeffs[i] * ip.base_polys[i](x) for i in eachindex(ip.coeffs))
end

function inter_poly_str(ip::InterpolationPolynomial)
    @argcheck length(ip.coeffs) == length(ip.base_polys) "interpolation polynomial has mismatched lengths"

    segments = String[]
    for i in 1:length(ip.coeffs)
        coeff = ip.coeffs[i]
        segment = ""
        if (coeff == 0)
            continue
        elseif (coeff == 1)
            segment = "[$(ip.base_polys[i].func_str)]"
        else
            segment = "$(ip.coeffs[i])[$(ip.base_polys[i].func_str)]"
        end
        push!(segments, segment)
    end
    
    return join(segments, " + ")
end

function create_lagrange(xs::Vector{<:Real}, ys::Vector{<:Real})
    @argcheck length(xs) == length(ys) "xs and ys must have equal lengths"
    @argcheck unique(xs) == xs "elements in xs must not appear more than once"
    @argcheck length(xs) >= 1 "there must be at least one point"

    coeffs = ys
    cardinal_funcs = _cardinal_functions(xs)

    return LagrangePolynomial(coeffs, cardinal_funcs)
end

function _cardinal_functions(xs::Vector{<:Real})
    len = length(xs)
    funcs = BasisPolynomial[]

    for n = 1:len
        # (x-xs[1])...(x-xs[n-1])(x-xs[n+1])...(x-xs[len]) evaluated at xs[n]
        denom = prod((xs[n] - xs[i]) for i in 1:len if i != n)
        # (x-xs[1])...(x-xs[n-1])(x-xs[n+1])...(x-xs[len])/(evaluated denom)
        func = (x -> prod((x - xs[i]) for i in 1:len if i != n) / denom)

        # build string of form [(x-xs[.])...(x-xs[.])]/(evaluated denom)
        segments = String[]
        absdenom = abs(denom)
        if (sign(denom) == -1)
            push!(segments, "-")
        end
        for i in 1:len
            if (i == n)
                continue
            elseif (sign(xs[i]) == 1)
                xsi = xs[i]
                segment = "(x - $xsi)"
                push!(segments, segment)
            elseif (sign(xs[i]) == -1)
                xsi = -xs[i]
                segment = "(x + $xsi)"
                push!(segments, segment)
            else
                segment = "x"
                push!(segments, segment)
            end
        end
        if (absdenom != 1)
            push!(segments, "/$absdenom")
        end
        string = join(segments)

        # instantiate a BasisPolynomial with corresponding function and string
        basis_poly = BasisPolynomial(func, string)
        push!(funcs, basis_poly)
    end

    return funcs
end

function create_newton(xs::Vector{<:Real}, ys::Vector{<:Real})
    @argcheck length(xs) == length(ys) "xs and ys must have equal lengths"
    @argcheck unique(xs) == xs "elements in xs must not appear more than once"
    @argcheck length(xs) >= 1 "there must be at least one point"
    
    len = length(xs)
    basis_funcs = _basis_functions(xs)
    coeffs = []
    push!(coeffs, ys[1])

    # compute coefficients
    for n = 2:len
        numer = ys[n] - sum(coeffs[n-i]*(basis_funcs[n-i](xs[n])) for i in 1:(n-1))
        denom = basis_funcs[n](xs[n])
        coeff = numer / denom
        push!(coeffs, coeff)
    end

    return NewtonPolynomial(coeffs, basis_funcs)
end

function _basis_functions(xs::Vector{<:Real})
    len = length(xs)
    funcs = BasisPolynomial[]
    
    # instantiate the first BasisPolynomial f1(x) = 1
    func = (_ -> 1)
    string = "1"
    basis_poly = BasisPolynomial(func, string)
    push!(funcs, basis_poly)

    for n = 2:len
        # recursively define the next function fn(x) = fn-1(x) * (x - xs[n-1])
        func = (x -> funcs[n-1](x) * (x - xs[n-1]))

        # build string of form "fn-1(x)(x - xs[n-1])"
        segments = String[]
        if (funcs[n-1].func_str != "1")
            push!(segments, funcs[n-1].func_str)
        end
        if (sign(xs[n-1]) == 1)
            segment = "(x - $(xs[n-1]))"
            push!(segments, segment)
        elseif (sign(xs[n-1]) == -1)
            segment = "(x + $(-xs[n-1]))"
            push!(segments, segment)
        else
            segment = "x"
            push!(segments, segment)
        end
        string = join(segments)

        # instantiate a BasisPolynomial with corresponding function and string
        basis_poly = BasisPolynomial(func, string)
        push!(funcs, basis_poly)
    end

    return funcs
end

function lagrange_to_standard()
    return
end

function newton_to_standard()
    return
end
