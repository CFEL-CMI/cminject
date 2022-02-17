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
    is_boundary_hit(boundary; (v=>getfield(blub, v) for v in fieldnames(typeof(phasePosition)))...)
end
