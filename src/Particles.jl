using LabelledArrays


## Abstract Particle definition

abstract type Particle end

function u(particle)
    particle._u
end

function p(particle)
    particle._p
end

# copied from LabelledArrays.jl, since it doesn't export this function
@inline _symnames(::Type{SLArray{S,T,N,L,Syms}}) where {S,T,N,L,Syms} = Syms
# helper for getting the type of the fields of a struct, from the struct type
@inline _get_fieldtype(t::Type{T}, f::Symbol) where T = t.types[findfirst(n->n==f, fieldnames(t))]

@generated function get_at(x::P,::Val{s}) where {s, P<:Particle}
    _utype = _get_fieldtype(x, :_u)
    _ptype = _get_fieldtype(x, :_p)
    if s ∈ _symnames(_utype)
        _eltype = eltype(_utype)
        quote
            Base.@_propagate_inbounds_meta
            Base.getproperty(x._u, s)::$_eltype
        end
    elseif s ∈ fieldnames(_ptype)
        _ftype = _get_fieldtype(_ptype, s)
        quote
            Base.@_propagate_inbounds_meta
            Base.getfield(x._p, s)::$_ftype
        end
    else
        quote
            Base.@_propagate_inbounds_meta
            Base.getfield(x, s)
        end
    end
end

Base.@propagate_inbounds function Base.getproperty(particle::P, s::Symbol) where P<:Particle
    get_at(particle,Val(s))
end


## Phase space position definitions

ParticlePhasePos2D{T} = @SLVector T (:x, :z, :vx, :vz)
ParticlePhasePos3D{T} = @SLVector T (:x, :y, :z, :vx, :vy, :vz)


## Abstract spherical particle definitions

abstract type SphericalParticle <: Particle end

function mass(particle::P) where P<:SphericalParticle
    particle.ρ * (4/3)particle.r^3 * π
end


## SphericalParticle2D definition

# TODO: auto-generate the following complex declarations via a readable macro DSL,
# e.g. a declaration like
#
# @particle struct SphericalParticle2D{T}
#    @phase_space T (:x, :z, :vx, :vz)
#    @properties (r::T, ρ::T)
#    @velocities (:x => :vx, :z => :vz)
# end

mutable struct SphericalParticle2DProps{T}
    r::T
    ρ::T
end

struct SphericalParticle2D{T} <: SphericalParticle
    _u::ParticlePhasePos2D{T}
    _p::SphericalParticle2DProps{T}

    SphericalParticle2D{T}(u, p) where T = let p = new{T}(u, p)
        # Force-realize the generated functions for access to all properties.
        # If we don't do this, using solve() will compile to slow versions, for whatever reason, and those will stick.
        # Note that these field accesses are compiled out (the compiler gets that they don't do anything during runtime)
        p._u; p._p; p.x; p.z; p.vx; p.vz; p.r; p.ρ;
        p
    end
    SphericalParticle2D{T}(; x, z, vx, vz, r, ρ) where T = SphericalParticle2D{T}(x, z, vx, vz, r, ρ)
    SphericalParticle2D{T}(x, z, vx, vz, r, ρ) where T = SphericalParticle2D{T}(
        ParticlePhasePos2D{T}(x, z, vx, vz),
        SphericalParticle2DProps{T}(r, ρ)
    )
end
# Defines a mapping of particle props type -> particle type
associated_particle_type(::Type{SphericalParticle2DProps{T}}) where T = SphericalParticle2D{T}


## Example instances

const example_u = ParticlePhasePos2D{Float64}(.1e-3, -0.1284, 7e-5, 3.63)
const example_p = SphericalParticle2DProps{Float64}(100e-9, 1050.0)
const example_particle = SphericalParticle2D{Float64}(example_u, example_p)