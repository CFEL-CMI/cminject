import PhysicalConstants.CODATA2018: k_B
const kB = k_B.val

include("Interpolation.jl")

# ------ Fields

abstract type Field
end

struct SumOfFields{N} <: Field
    fields::NTuple{N, Field}
end

function Base.:+(f1::Field, f2::Field)
    SumOfFields{2}((f1, f2))
end

function acceleration(particle::Particle, field::SumOfFields, time)
    a = zeros(2)  # FIXME
    for subfield in field.fields
        a .+= acceleration(particle, subfield, time)
    end
end

struct StokesFlowField{T,D,A} <: Field
    # TODO replace with struct from Interpolations.jl
    grid::RegularGrid{D,3,T,A}
    T::T
    mᶠ::T
    μ::T
end

const example_grid = hdf5_to_regulargrid(
    "/home/chip/Documents/cminject-materials/goldspheres-27nm/nitrogen_1.8mbar_extended_nan.h5"
)
const example_field = StokesFlowField(example_grid, 293.15, 4.27e-26, 1.76e-5)

function stokes_vars(
    particle::SphericalParticle{Float64},
    field::StokesFlowField{Float64, 2},
)
    u, P = particle._u, particle._params
    rₚ, mₚ = P.r, mass(particle)

    vˣᶠ, vᶻᶠ, p = interpol1(field.grid, u.x, u.z)
    Knf = field.μ * √(π * kB * field.T / (2*field.mᶠ))
    Kn = Knf / (p * rₚ)
    Cc = 1 + Kn*(1.231 + 0.4695exp(-1.1783 / Kn))
    a₀ = 6π * field.μ * P.r / Cc / mₚ
    Δvˣ, Δvᶻ = (vˣᶠ - u.vx), (vᶻᶠ - u.vz)

    return isnan(p), a₀, Δvˣ, Δvᶻ
end

function acceleration(
    particle::SphericalParticle{Float64},
    field::StokesFlowField{Float64, 2},
    time
)
    invalid, a₀, Δvˣ, Δvᶻ = stokes_vars(particle, field)
    invalid ? (NaN, NaN) : (a₀ * Δvˣ, a₀ * Δvᶻ)
end