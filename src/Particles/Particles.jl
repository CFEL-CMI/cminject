using LabelledArrays

include("ParticlesBase.jl")


## SphericalParticle2D definition
@declare_particle struct SphericalParticle2D{T}
    @phase_space T (:x, :z, :vx, :vz)
    @velocities (:x => :vx, :z => :vz)
    r::T
    ρ::T
end


## Examples
const example_particle = SphericalParticle2D{Float64}(r=100e-9, ρ=1050.0, x=1e-4, z=-0.1284, vx=7e-5, vz=3.63)
const example_u = example_particle._u
const example_p = example_particle._p