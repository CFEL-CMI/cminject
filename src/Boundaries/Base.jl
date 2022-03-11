"""
Abstract supertype for all boundary types in CMInject.

Subtypes `MyBoundary <: AbstractBoundary` should at least contain the following fields:
- ε: The tolerance

They should also implement the following methods:
- is_boundary_hit(boundary::MyBoundary)
"""
abstract type AbstractBoundary end

"""
    is_boundary_hit(boundary, phasePosition)

True if the `phasePosition` is considered to have hit the `boundary`.
Converts the `phasePosition` labelled array into function arguments for the internal methods.
"""
function is_boundary_hit(boundary::B, phasePosition::LArray) where B<:AbstractBoundary
    if (size(phasePosition)[1] == 6)
        is_boundary_hit(boundary; x=phasePosition.x, y=phasePosition.y, z=phasePosition.z)
    elseif (size(phasePosition)[1] == 4)
        is_boundary_hit(boundary; y=phasePosition.y, z=phasePosition.z)
    else
        error("Dimensions $(size(phasePosition)) not (yet) supported.")
    end
end
