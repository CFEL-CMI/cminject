"""
    Skimmer

Represents a skimmer at position (x),y,z with radius r
"""
struct Skimmer{T} <: AbstractBoundary where T<:Number
    x::T
    y::T
    z::T
    r::T
    ε::T
end

function is_boundary_hit(boundary::Skimmer; x=boundary.x, y=boundary.y, z=boundary.z)
    return z ≥ boundary.z-boundary.ε && z ≤ boundary.z &&
        √((x-boundary.x)^2 + (y-boundary.y)^2) ≥ boundary.r
end

