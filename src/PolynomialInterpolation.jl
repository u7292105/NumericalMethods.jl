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
    return sum(ip.coeffs[i] * p.base_polys[i](x) for i in eachindex(p.coeffs))
end

function inter_poly_str(ip::InterpolationPolynomial)
    return ""
end

function create_lagrange()
    return
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

