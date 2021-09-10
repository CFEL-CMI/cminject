#using Lazy
#using StaticArrays
using LabelledArrays

# ----- Particles

abstract type Particle{T}
end

function get_u(p::P where P<:Particle)
    p._u
end

function get_params(p::P where P<:Particle)
    p._params
end

function Base.setproperty!(particle::P where P<:Particle, f::Symbol, v)
    # provide a less confusing error than "particle does not have field x" or such
    error(
        "Particles are immutable, cannot set a property! Construct a new Particle instead."
    )
end

struct _SphericalParticleProps{T}
    ρ::T
    r::T
    affe::Int64
end

struct SphericalParticle{T} <: Particle{T}
    _params::_SphericalParticleProps{T}
    _u::@SLVector T (:x, :z, :vx, :vz)
end
SphericalParticle(ρ::T, r::T, x::T, z::T, vx::T, vz::T) where T = (
    SphericalParticle{T}(_SphericalParticleProps(ρ, r, 12), SLVector(x=x, z=z, vx=vx, vz=vz))
)

const example_particle = SphericalParticle(1050.0, 100e-9, 0.0, -0.128, 0.0, 10.0)

Base.@propagate_inbounds function Base.getproperty(particle::SphericalParticle{T}, f::Symbol) where T
    if f ∈ (:x, :z, :vx, :vz)
        Base.getproperty(particle._u, f)
    elseif f ∈ (:ρ, :r, :affe)
        Base.getfield(particle._params, f)
    else
        Base.getfield(particle, f)  # required so `particle._u` and `particle._params` can resolve
    end
end

function u(particle::SphericalParticle)
    particle.u
end

function mass(particle::SphericalParticle)
    particle._params.ρ * (4/3)particle._params.r^3 * π
end