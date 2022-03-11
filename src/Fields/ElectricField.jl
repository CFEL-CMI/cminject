# Note that you can't directly include this file, you'll have to include "Fields.jl"

"""
    ElectricField2d(; norm, gradients)

Representation of a 2d electric field.
The data is stored in an interpolated form,
not merely the values at the grid points are available.

It consists of:
- `norm`, the (two dimensional) interpolator for the norm of the field
- `gradients`, either the two gradient interpolators or nothing,
    if the gradients should be calculated via the `norm`
"""
struct ElectricField2D{ITP<:AbstractInterpolation} <: Field
    norm::ITP
    gradients::Union{Some{Tuple{ITP,ITP}}, Nothing}
end
function ElectricField2D(itp::ITP, gradients::Tuple{ITP,ITP}) where ITP
    ElectricField2D{ITP}(itp, Some(gradients))
end
function ElectricField2D(itp::ITP) where ITP
    ElectricField2D{ITP}(itp, nothing)
end
function ElectricField2D(hdf5Path::s) where s<:AbstractString
    ElectricField2D(interpolateElectricField(hdf5Path)...)
end
Base.show(io::IO, f::ElectricField2D) = print(io, "ElectricField2D")

"""
    ElectricField(; norm, gradients)

Representation of an electric field.
The data is stored in an interpolated form,
not merely the values at the grid points are available.

It consists of:
- `norm`, the (three dimensional) interpolator for the norm of the field
- `gradients`, either the three gradient interpolators or nothing,
    if the gradients should be calculated via the `norm`
"""
struct ElectricField{ITP<:AbstractInterpolation} <: Field
    norm::ITP
    gradients::Union{Some{Tuple{ITP,ITP,ITP}}, Nothing}
end
function ElectricField(itp::ITP, gradients::Tuple{ITP,ITP,ITP}) where ITP
    ElectricField{ITP}(itp, Some(gradients))
end
function ElectricField(itp::ITP) where ITP
    ElectricField{ITP}(itp, nothing)
end
function ElectricField(hdf5Path::s) where s<:AbstractString
    ElectricField(interpolateElectricField(hdf5Path)...)
end
Base.show(io::IO, f::ElectricField) = print(io, "ElectricField")

function acceleration(particle, field::ElectricField2D, time)
    force = -getEnergyGradient(field, particle.starkCurve, particle.x, particle.y)
    acceleration = force/mass(particle)
    (vx = acceleration[1], vy = acceleration[2])
end

function acceleration(particle, field::ElectricField, time)
    # TODO: Do something such that the number of dimensions is not pre-determined
    force = -getEnergyGradient(field, particle.starkCurve, particle.x, particle.y, particle.z)
    acceleration = force/mass(particle)
    (vx = acceleration[1], vy = acceleration[2], vz = acceleration[3])
end

"""
    getEnergyGradient(field, starkCurve, x, y, z)

Calculates the energy gradient of the specified electric field
for a particle with the specified stark curve at the specified point in space.

The arguments are:
- `field`, the electric field the particle is in
- `starkCurve`, the stark curve for the particle
- `x,y,z`, the point in space of the particle
"""
function getEnergyGradient(field, starkCurve, x::T, y::T, z::T)::Vector{T} where T<:Number
    # The gradient of the energy over space is given by
    # ∇E(r) = E'(ε(r))⋅∇ε(r)
    if (field.gradients === nothing)
        gradient(starkCurve, field.norm(x,y,z)) .* gradient(field.norm, x,y,z)
    else
        gradients = [gradient(x,y,z) for gradient ∈ field.gradients.value]
        gradient(starkCurve, field.norm(x,y,z)) .* gradients
    end
end

"""
    getEnergyGradient(field, starkCurve, x, y)

Calculates the energy gradient of the specified electric field
for a particle with the specified stark curve at the specified point in space.

The arguments are:
- `field`, the electric field the particle is in
- `starkCurve`, the stark curve for the particle
- `x,y`, the point in space of the particle
"""
function getEnergyGradient(field, starkCurve, x::T, y::T)::Vector{T} where T<:Number
    # The gradient of the energy over space is given by
    # ∇E(r) = E'(ε(r))⋅∇ε(r)
    if (field.gradients === nothing)
        gradient(starkCurve, field.norm(x,y)) .* gradient(field.norm, x,y)
    else
        gradients = [gradient(x,y) for gradient ∈ field.gradients.value]
        gradient(starkCurve, field.norm(x,y)) .* gradients
    end
end
