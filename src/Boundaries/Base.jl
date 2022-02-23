"""
Abstract supertype for all boundary types in CMInject.

Subtypes `MyBoundary <: AbstractBoundary` should at least contain the following fields:
- ε: The tolerance

They should also implement the following methods:
- is_boundary_hit(boundary::MyBoundary)
"""
abstract type AbstractBoundary end

function is_boundary_hit(boundary::B, phasePosition::LArray) where B<:AbstractBoundary
    asTuple = convert(NamedTuple, phasePosition)
    # Heavily inspired by https://discourse.julialang.org/t/get-fieldnames-and-values-of-struct-as-namedtuple/8991
    is_boundary_hit(boundary; (v=>getfield(asTuple, v) 
                               for v ∈ fieldnames(typeof(asTuple)) 
                               # Filter out velocities
                               if String(v)[1] != 'v')...)
end
