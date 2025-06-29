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
    return sum(ip.coeffs[i] * ip.base_polys[i](x) for i in eachindex(p.coeffs))
end

function inter_poly_str(ip::InterpolationPolynomial)
    return ""
end

function create_lagrange()
    return
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
        if (sign(denom) == -1)
            push!(segments, "-")
            denom = -denom
        end
        push!(segments, "[")
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
        push!(segments, "]")
        if (denom > 1)
            push!(segments, "/$denom")
        end
        string = join(segments)

        # instantiate a BasisPolynomial with corresponding function and string
        basis_poly = BasisPolynomial(func, string)
        push!(funcs, basis_poly)
    end

    return funcs
end

function create_newton()
    return
end

function lagrange_to_standard()
    return
end

function newton_to_standard()
    return
end
