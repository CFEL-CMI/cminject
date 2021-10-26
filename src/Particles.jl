using LabelledArrays

include("ParticlesBase.jl")
include("Interpolation.jl")


## SphericalParticle2D definition
@declare_particle struct SphericalParticle2D{T}
    @phase_space T (:x, :z, :vx, :vz)
    @velocities (:x => :vx, :z => :vz)
    r::T
    ρ::T
end

#=
    StarkParticle2d

(Neutral) 2d particle that interacts with electrical fields via the Stark effect.
The stark curve is filled in upon creation from the source.
Hence the source also decides which _specific_ molecule this is.
=#
# TODO: Restrict ITP to AbstractInterpolation
@declare_particle struct StarkParticle2d{T, ITP}
    @phase_space T (:x, :y, :vx, :vy)
    @velocities (:x => :vx, :y => :vy)
    # Mass
    m::T
    starkCurve::ITP
end

#=
    StarkParticle

(Neutral) particle that interacts with electrical fields via the Stark effect.
The stark curve is filled in upon creation from the source.
Hence the source also decides which _specific_ molecule this is.
=#
# TODO: Code duplicate
@declare_particle struct StarkParticle{T, ITP}
    @phase_space T (:x, :y, :z, :vx, :vy, :vz)
    @velocities (:x => :vx, :y => :vy, :z => :vz)
    # Mass
    m::T
    starkCurve::ITP
end

function mass(particle::StarkParticle2d)
    particle.m
end

function mass(particle::StarkParticle)
    particle.m
end

## Examples
const example_particle = SphericalParticle2D{Float64}(r=100e-9, ρ=1050.0, x=1e-4, z=-0.1284, vx=7e-5, vz=3.63)
const example_u = example_particle._u
const example_p = example_particle._p
