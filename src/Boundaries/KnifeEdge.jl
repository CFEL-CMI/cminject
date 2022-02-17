"""
    KnifeEdge

Represents a knife edge at position y,z in positive (p=1) or negative (p=-1) direction
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
