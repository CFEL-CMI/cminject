"""
Abstract supertype for all boundary types in CMInject.

Subtypes `MyBoundary <: AbstractBoundary` should at least contain the following fields:
- ε: The tolerance

They should also implement the following methods:
- is_boundary_hit(boundary::MyBoundary)
"""
abstract type AbstractBoundary end

# TODO: Specify type of phaseposition
function is_boundary_hit(boundary::B, phasePosition) where B<:AbstractBoundary
    # Heavily inspired by https://discourse.julialang.org/t/get-fieldnames-and-values-of-struct-as-namedtuple/8991
    call_hit_function(boundary; (v=>getfield(phasePosition, v) for v in fieldnames(typeof(phasePosition)))...)
end

function call_hit_function(boundary::B; __x::Vector{T}=Float64[]) where {B<:AbstractBoundary, T<:Number}
    if (size(__x)[1] == 6)
        # Three dimensions
        is_boundary_hit(boundary; x=__x[1], y=__x[2], z=__x[3])
    elseif (size(__x)[1] == 4)
        # Four dimensions
        is_boundary_hit(boundary; y=__x[1], z=__x[2])
    else
        error("Dimensions: $(size(__x)[1]) not supported")
    end
end
