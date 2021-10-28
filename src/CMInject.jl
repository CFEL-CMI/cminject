module CMInject

include("Interpolation.jl")
include("Particles.jl")
include("Fields.jl")
include("Sources.jl")
include("Detectors.jl")
include("ResultPlots.jl")
include("Experiment.jl")

using DifferentialEquations

## Benchmarking utilities - not to actually solve physical problems

const _du = similar(CMInject.example_particle._u)
function _quick()
    p = CMInject.example_particle
    f = CMInject.example_field
    f!(_du, p._u, (p._p, [p._p,], (f,)), 0.0)
    _du
end


## Example experiment

const example_dists = (
    x = Normal(0.0, 1e-3), z = Dirac(-0.128),
    vx = Normal(0.0, 0.1), vz = Dirac(10.0),
    r = Normal(13.5e-9, 2e-9), ρ = Dirac(1050.0)
)
const example_source = SamplingSource{SphericalParticle2D{Float64}}(example_dists)
const example_experiment = CMInject.Experiment(;
    source=CMInject.example_source,
    n_particles=1000,
    fields=(CMInject.example_field,),
    detectors=Tuple(CMInject.SectionDetector{Float64,:z}.(-0.003:0.001:0.003, true)),
    solver=EulerHeun(), time_span=(0.0, 0.03), time_step=1e-5, ensemble_alg=EnsembleThreads()
)

end