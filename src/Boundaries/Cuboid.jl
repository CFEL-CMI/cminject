"""
    Cuboid(; x1, x2, y1, y2, z1, z2, ε)

Represents a cuboid (e.g. a detector) from `x1`-`x2`, `y1`-`y2` and `z1`-`z2`
with open ends at `z1` and `z2` with tolerance `ε`.
"""
struct Cuboid{T} <: AbstractBoundary where T<:Number
    x1::T
    x2::T
    y1::T
    y2::T
    z1::T
    z2::T
    ε::T
end

function is_boundary_hit(boundary::Cuboid;
        x=(boundary.x1+boundary.x2)/2,
        y=(boundary.y1+boundary.y2)/2,
        z=(boundary.z1+boundary.z2)/2)
    return z ≥ boundary.z1 && z ≤ boundary.z2 &&
        (x ≤ boundary.x1 || x ≥ boundary.x2 ||
         y ≤ boundary.y1 || y ≥ boundary.y2)
end

