import PhysicalConstants.CODATA2018: k_B
include("StarkEffect.jl")
using .StarkEffect
const kB = k_B.val

include("Interpolation.jl")

# ------ Fields

abstract type Field
end
noise(particle, f::Field, time) = ()  # no noise term for fields required by default


struct StokesFlowField{T,ITP<:AbstractInterpolation} <: Field
    interpolator::ITP
    T::T
    mᶠ::T
    μ::T
end
function StokesFlowField(itp::ITP, t::T, mf::T, mu::T) where {ITP, T}
    StokesFlowField{T,ITP}(itp, t, mf, mu)
end
Base.show(io::IO, f::StokesFlowField) = print(io, "StokesFlowField(T=$(f.T),mᶠ=$(f.mᶠ),μ=$(f.μ))")

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

const example_itp = hdf5_to_interpolator(joinpath(@__DIR__, "../data/nitrogen_1.8mbar_extended_nan.h5"))
const example_field = StokesFlowField(example_itp, 293.15, 4.27e-26, 1.76e-5)

function stokes_vars(particle, field)
    rₚ, mₚ = particle.r, mass(particle)
    xₚ, zₚ, vˣₚ, vᶻₚ = particle.x, particle.z, particle.vx, particle.vz
    μ, T, mᶠ = field.μ, field.T, field.mᶠ

    vˣᶠ, vᶻᶠ, p = field.interpolator(xₚ, zₚ)
    Kn = μ * √(π * kB * T / 2mᶠ) / (p * rₚ)
    # TODO this Cunningham correction is only one specific model. Make it parameterizable (-> field property)
    Cc = 1 + Kn*(1.231 + 0.4695exp(-1.1783 / Kn))
    a₀ = 6π * μ * rₚ / (Cc * mₚ)
    Δvˣ, Δvᶻ = (vˣᶠ - vˣₚ), (vᶻᶠ - vᶻₚ)

    return isnan(p), a₀, Δvˣ, Δvᶻ, Cc
end

function acceleration(particle, field::StokesFlowField, time)
    invalid, a₀, Δvˣ, Δvᶻ, _ = stokes_vars(particle, field)
    invalid ? (vx=NaN, vz=NaN) : (vx = a₀ * Δvˣ, vz = a₀ * Δvᶻ)
end

function noise(particle, field::StokesFlowField, time)
    invalid, _, _, _, Cc = stokes_vars(particle, field)
    s₀fπ = π * 216field.μ * kB * field.T / (π^2 * (2particle.r)^5 * particle.ρ^2)
    s₀ = √(s₀fπ / Cc)
    invalid ? (vx=0.0, vz=0.0) : (vx=s₀, vz=s₀)
end

function acceleration(particle, field::ElectricField, time)
    # TODO: Obtain Stark curve from particle
    starkCurve=calculateStarkCurves(1e-3, 0.0, 1e-2, 0, 5, 0, 0, 9.93910344e-26, 3.31303448e-26, 1, 1, 1, 1.66782048e-29)[5]
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
