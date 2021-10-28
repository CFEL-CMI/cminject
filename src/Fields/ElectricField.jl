# Note that you can't directly include this file, you'll have to include "Fields.jl"
include("../Physics/StarkEffect.jl")
using .StarkEffect

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

function acceleration(particle, field::ElectricField, time)
    # TODO: Obtain Stark curve from particle
    allStarkCurves = calculateStarkCurves(1e-3, 0.0, 1e-2, 0, 5, 0, 0, 9.93910344e-26, 3.31303448e-26, 1, 1, 1, 1.66782048e-29)
    @inbounds starkCurve=allStarkCurves[5]
    # TODO: Obtain location from particle
    r = [1.5, 1.2]
    force = getEnergyGradient(field.interpolator, starkCurve, r)
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
