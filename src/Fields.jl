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

#function acceleration(particle, field::SumOfFields, time)
#    a = zeros(2)  # FIXME
#    for subfield in field.fields
#        a .+= acceleration(particle, subfield, time)
#    end
#end

struct StokesFlowField{T,ITP<:AbstractInterpolation} <: Field
    # TODO replace with struct from Interpolations.jl
    interpolator::ITP
    T::T
    mᶠ::T
    μ::T
end
function StokesFlowField(itp::ITP, t::T, mf::T, mu::T) where {ITP, T}
    StokesFlowField{T,ITP}(itp, t, mf, mu)
end
Base.show(io::IO, f::StokesFlowField) = print(io, "StokesFlowField(T=$(f.T),mᶠ=$(f.mᶠ),μ=$(f.μ))")

const example_itp = hdf5_to_interpolator(
    "/home/chip/Documents/cminject-materials/goldspheres-27nm/nitrogen_1.8mbar_extended_nan.h5"
)
const example_field = StokesFlowField(example_itp, 293.15, 4.27e-26, 1.76e-5)

function stokes_vars(particle, field)
    rₚ, mₚ = particle.r, particle.ρ * 0.75particle.r^3 * π

    vˣᶠ, vᶻᶠ, p = field.interpolator(particle.x, particle.z)
    Knf = field.μ * √(π * kB * field.T / (2*field.mᶠ))
    Kn = Knf / (p * rₚ)
    Cc = 1 + Kn*(1.231 + 0.4695exp(-1.1783 / Kn))
    a₀ = 6π * field.μ * rₚ / Cc / mₚ
    Δvˣ, Δvᶻ = (vˣᶠ - particle.vx), (vᶻᶠ - particle.vz)

    return isnan(p), a₀, Δvˣ, Δvᶻ
end

function stokes_vars_direct(particle, field)
    rₚ, mₚ = particle._p.r, particle._p.ρ * 0.75particle._p.r^3 * π

    vˣᶠ, vᶻᶠ, p = field.interpolator(particle._u.x, particle._u.z)
    Knf = field.μ * √(π * kB * field.T / (2*field.mᶠ))
    Kn = Knf / (p * rₚ)
    Cc = 1 + Kn*(1.231 + 0.4695exp(-1.1783 / Kn))
    a₀ = 6π * field.μ * rₚ / Cc / mₚ
    Δvˣ, Δvᶻ = (vˣᶠ - particle._u.vx), (vᶻᶠ - particle._u.vz)

    return isnan(p), a₀, Δvˣ, Δvᶻ
end

function acceleration(
    particle,
    field,
    time
) where T
    invalid, a₀, Δvˣ, Δvᶻ = stokes_vars(particle, field)
    invalid ? (NaN, NaN) : (a₀ * Δvˣ, a₀ * Δvᶻ)
end

function acceleration_direct(
    particle,
    field,
    time
) where T
    invalid, a₀, Δvˣ, Δvᶻ = stokes_vars_direct(particle, field)
    invalid ? (NaN, NaN) : (a₀ * Δvˣ, a₀ * Δvᶻ)
end