#=
    StarkParticle2D

(Neutral) 2d particle that interacts with electrical fields via the Stark effect.
The stark curve is filled in upon creation from the source.
Hence the source also decides which _specific_ molecule this is.
=#
# TODO: Restrict ITP to AbstractInterpolation
@declare_particle struct StarkParticle2D{T,ITP}
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
@declare_particle struct StarkParticle{T,ITP}
    @phase_space T (:x, :y, :z, :vx, :vy, :vz)
    @velocities (:x => :vx, :y => :vy, :z => :vz)
    # Mass
    m::T
    starkCurve::ITP
end

function mass(particle::Union{StarkParticle2D,StarkParticle})
    particle.m
end
