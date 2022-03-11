"""
    KnifeEdge(; y, z, p, ε)

Represents a knife edge at position y,z in positive (p=1) or negative (p=-1) direction

It consists of:
- `y`, the y position (*height*) of the knife edge
- `z`, the z position of the knife edge
- `p`, whether it extends to the negative (`p<0`) or positive (`p>0`) y direction
- `ε`, the tolerance for `z`
"""
struct KnifeEdge{T,B} <: AbstractBoundary where {T<:Number, B<:Signed}
    y::T
    z::T
    p::B
    ε::T
end

function is_boundary_hit(boundary::KnifeEdge; x=0, y=boundary.y, z=boundary.z)
    return z ≥ boundary.z-boundary.ε && z ≤ boundary.z &&
        (y ≤ boundary.y && boundary.p < 0 || y ≥ boundary.y && boundary.p > 0)
end
