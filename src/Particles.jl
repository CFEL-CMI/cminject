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
end

struct SphericalParticle{T} <: Particle{T}
    _params::_SphericalParticleProps{T}
    _u::@SLVector T (:x, :z, :vx, :vz)
end
SphericalParticle(ρ::T, r::T, x::T, z::T, vx::T, vz::T) where T = (
    SphericalParticle{T}(_SphericalParticleProps(ρ, r), SLVector(x=x, z=z, vx=vx, vz=vz))
)
#SphericalParticle(ρ, r, affe::Int32, x, z, vx, vz) = let x_, z_, vx_, vz_ = promote(x, z, vx, vz)
#    SphericalParticle(_SphericalParticleProps(ρ, r, affe), SLVector(x=x_, z=z_, vx=vx_, vz=vz_))
#end

const example_particle = SphericalParticle(1050.0, 100e-9, 0.0, -0.128, 0.0, 10.0)


#= seems to be inefficient... better do via some compiler/macro tricks?
Base.@propagate_inbounds function Base.getproperty(particle::SphericalParticle, f::Symbol)
    if f in (:x, :z, :vx, :vz)
        Base.getproperty(particle._u, f)
    elseif f in (:ρ, :r, :affe)
        Base.getfield(particle._p, f)
    else
        Base.getfield(particle, f)
    end
end
=#



#=
function u(particle::P where P <: Particle)
    particle.u
end

struct SphericalParticle{T<:Real, N} <: Particle{T}
    u::Vector{T}
    r::T
    ρ::T
    T::T
end

SphericalParticle(u::Vector{T}, r::T, ρ::T) where T<:Real = SphericalParticle{eltype(u), length(u)}(u, r, ρ)

@forward SphericalParticle.u (
    Base.:+, Base.:*, Base.:/,
    Base.similar, Base.size, Base.length,
    Base.iterate, Base.getindex, Base.setindex!
)
=#


function u(particle::SphericalParticle)
    particle.u
end

function mass(particle::SphericalParticle)
    particle._params.ρ * (4/3)particle._params.r^3 * π
end