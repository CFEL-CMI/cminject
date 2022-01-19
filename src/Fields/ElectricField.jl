# Note that you can't directly include this file, you'll have to include "Fields.jl"

"""
    ElectricField2d

Representation of a 2d electric field.
The data is stored in an interpolated form,
not merely the values at the grid points are available.
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
Base.show(io::IO, f::ElectricField2D) = print(io, "ElectricField2D")

"""
    ElectricField

Representation of an electric field.
The data is stored in an interpolated form,
not merely the values at the grid points are available.
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
Base.show(io::IO, f::ElectricField) = print(io, "ElectricField")

function acceleration(particle, field::ElectricField2D, time)
    r = [particle.x, particle.y]
    force = getEnergyGradient(field, particle.starkCurve, r)
    acceleration = force/mass(particle)
    (vx = acceleration[1], vy = acceleration[2])
end

function acceleration(particle, field::ElectricField, time)
    # TODO: Do something such that the number of dimensions is not pre-determined
    r = [particle.x, particle.y, particle.z]
    force = getEnergyGradient(field, particle.starkCurve, r)
    acceleration = force/mass(particle)
    (vx = acceleration[1], vy = acceleration[2], vz = acceleration[3])
end

"""
    getEnergyGradient(field, starkCurve, r)

Calculates the energy gradient of the specified electric field
"""
function getEnergyGradient(field, starkCurve, r)::AbstractArray
    # The gradient of the energy over space is given by
    # ∇E(r) = E'(ε(r))⋅∇ε(r)
    if (field.gradients == nothing)
        gradient(starkCurve, field.norm(r...)) .* gradient(field.norm, r...)
    else
        gradient(starkCurve, field.norm(r...)) .* [gradient(r...) for gradient ∈ field.gradients.value]
    end
end
