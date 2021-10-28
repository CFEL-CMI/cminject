#=
    SphericalParticle2D

A spherical particle moving in two dimensions `x` and `z` with associated velocities `vx` and `vz`,
with an associated radius `r` and a material density `ρ`.
=#
@declare_particle struct SphericalParticle2D{T}
    @phase_space T (:x, :z, :vx, :vz)
    @velocities (:x => :vx, :z => :vz)
    r::T
    ρ::T
end


#=
    SphericalParticle3D

A spherical particle moving in three dimensions `x`, `y`, `z` with associated velocities `vx`, `vy`, `vz`,
with an associated radius `r` and a material density `ρ`.
=#
@declare_particle struct SphericalParticle3D{T}
    @phase_space T (:x, :y, :z, :vx, :vy, :vz)
    @velocities (:x => :vx, :y => :vy, :z => :vz)
    r::T
    ρ::T
end