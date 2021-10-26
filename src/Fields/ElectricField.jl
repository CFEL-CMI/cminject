# Note that you can't directly include this file, you'll have to include "Fields.jl"
include("../Physics/StarkEffect.jl")
using .StarkEffect
include("../Interpolation.jl")

"""
    ElectricField2d

Representation of a 2d electric field.
The data is stored in an interpolated form,
not merely the values at the grid points are available.
"""
struct ElectricField2d{ITP<:AbstractInterpolation} <: Field
    interpolator::ITP
end
function ElectricField2d(itp::ITP) where ITP
    ElectricField{ITP}(itp)
end
Base.show(io::IO, f::ElectricField2d) = print(io, "ElectricField2d")

"""
    ElectricField

Representation of an electric field.
The data is stored in an interpolated form,
not merely the values at the grid points are available.
"""
struct ElectricField{ITP<:AbstractInterpolation} <: Field
    interpolator::ITP
end
function ElectricField(itp::ITP) where ITP
    ElectricField{ITP}(itp)
end
Base.show(io::IO, f::ElectricField) = print(io, "ElectricField")

function acceleration(particle, field::ElectricField2d, time)
    r = [particle.x, particle.y]
    force = getEnergyGradient(field.interpolator, particle.starkCurve, r)
    force/mass(particle)
end

function acceleration(particle, field::ElectricField, time)
    # TODO: Do something such that the number of dimensions is not pre-determined
    r = [particle.x, particle.y, particle.z]
    force = getEnergyGradient(field.interpolator, particle.starkCurve, r)
    force/mass(particle)
end

"""
    getEnergyGradient(field, starkCurve, r)

Calculates the energy gradient of the specified electric field
"""
function getEnergyGradient(field, starkCurve, r)::AbstractArray
    # The gradient of the energy over space is given by
    # ∇E(r) = E'(ε(r))⋅∇ε(r)
    gradient(starkCurve, field(r...)) .* gradient(field, r...)
end
