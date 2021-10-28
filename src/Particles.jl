using LabelledArrays

include("Particles/Base.jl")
include("Particles/ClassicalParticles.jl")
include("Particles/StarkParticles.jl")


## Examples
const example_particle = SphericalParticle2D{Float64}(r=100e-9, ρ=1050.0, x=1e-4, z=-0.1284, vx=7e-5, vz=3.63)
const example_u = example_particle._u
const example_p = example_particle._p
