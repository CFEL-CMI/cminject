#using Lazy
#using StaticArrays
using LabelledArrays

# ----- Particles

#abstract type Particle{T} where T end

mutable struct SphericalParticleProps{T}
    r::T
    ρ::T
end
ParticlePhasePos2D{T} = @SLVector T (:x, :z, :vx, :vz)
struct SphericalParticle{T}
    _u::ParticlePhasePos2D{T}
    _p::SphericalParticleProps{T}
end
SphericalParticle{T}(x, z, vx, vz, r, ρ) where T = SphericalParticle{T}(
    ParticlePhasePos2D{T}(x, z, vx, vz),
    SphericalParticleProps{T}(r, ρ)
)

@inline @generated function Base.getindex(x::SphericalParticle{T},::Val{s}) where {s, T}
    if s ∈ (:x, :z, :vx, :vz)
        :(x._u[s]::T)
    elseif s ∈ (:ρ, :r)
        _ptype = x.types[findfirst(n->n==:_p, fieldnames(x))]
        _ftype = _ptype.types[findfirst(n->n==s, fieldnames(_ptype))]
        :(x._p.$s::$_ftype)
    else
        :(Base.getfield(x, s))
    end
end

@inline function Base.getproperty(particle::SphericalParticle{T}, s::Symbol) where T
    getindex(particle,Val(s))
    #= if s ∈ (:x, :z, :vx, :vz)
        :(Base.getproperty(particle._u, s))
    elseif s ∈ (:ρ, :r)
        :(Base.getfield(particle._p, s))
    else
        :(Base.getfield(particle, s))
    end =#
end

const example_u = ParticlePhasePos2D{Float64}(.1e-3, -0.1284, 7e-5, 3.63)
const example_p = SphericalParticleProps{Float64}(100e-9, 1050.0)
const example_particle = SphericalParticle{Float64}(example_u, example_p)


function u(particle)
    particle.u
end

function mass(particle)
    particle.p.ρ * (4/3)particle.p.r^3 * π
end