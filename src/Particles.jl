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

@inline @generated function get_at(x::SphericalParticle{T},::Val{s}) where {s, T}
    _utype = x.types[findfirst(n->n==:_u, fieldnames(x))]
    _ptype = x.types[findfirst(n->n==:_p, fieldnames(x))]
    if s ∈ _utype.parameters[5]
        _eltype = _utype.parameters[2]
        :(Base.getproperty(x._u, s)::$_eltype)
    elseif s ∈ fieldnames(_ptype)
        _ftype = _ptype.types[findfirst(n->n==s, fieldnames(_ptype))]
        :(Base.getfield(x._p, s)::$_ftype)
    else
        :(Base.getfield(x, s))
    end
end

@inline function Base.getproperty(particle::SphericalParticle{T}, s::Symbol) where T
    get_at(particle,Val(s))
end

const example_u = ParticlePhasePos2D{Float64}(.1e-3, -0.1284, 7e-5, 3.63)
const example_p = SphericalParticleProps{Float64}(100e-9, 1050.0)
const example_particle = SphericalParticle{Float64}(example_u, example_p)
let p = example_particle
    # force-realize the generated functions for access to all properties
    # if we don't do this, using @btime will compile to slow versions,
    # for whatever reason, and those will stick...
    p.x; p.z; p.vx; p.vz; p.r; p.ρ;
end;


function u(particle)
    particle.u
end

function mass(particle)
    particle.p.ρ * (4/3)particle.p.r^3 * π
end