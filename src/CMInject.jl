module CMInject

export Experiment, simulate,
    SamplingSource,
    StokesFlowField,
    ElectricField, ElectricField2D,
    SectionDetector,
    SphericalParticle2D, SphericalParticle3D,
    StarkParticle, StarkParticle2D,
    HDF5ResultStorage, store_results!, plot_results

include("Interpolation.jl")
include("Particles.jl")
include("Fields.jl")
include("Sources.jl")
include("Detectors.jl")
include("ResultPlots.jl")
include("Experiment.jl")
include("ResultStorage.jl")

using DifferentialEquations

## Benchmarking utilities - not to actually solve physical problems

const _du = similar(CMInject.example_particle._u)
function _quick()
    p = CMInject.example_particle
    f = CMInject.example_field
    f!(_du, p._u, (p._p, [p._p,], (f,)), 0.0)
    _du
end


## Example instances to play around with

"""
    get_example_instances()

Returns a NamedTuple containing example instances of some CMInject, including a full working Experiment:

- `field`: A StokesFlowField of an ALS
- `field_interpolator`: An interpolator (see Interpolations.jl) that's used by `field`
- `source`: A `SamplingSource` generating `SphericalParticle2D{Float64}` instances, based on `source_distributions`
- `source_distributions`: A set of distributions for all properties making up `SphericalParticle2D{Float64}` instances
- `detectors`: A Tuple of `SectionDetector{Float64, :z}`s at z positions (-3e-3, -2e-3, -1e-3, 0, 1e-3, 2e-3, 3e-3)
- `experiment`: An Experiment tying all of the above together and simulating 1000 particles for the timespan
  t = [0s, 0.05s].
"""
function get_example_instances()
    example_flow_field_file = joinpath(@__DIR__, "../data/nitrogen_1.8mbar_extended_nan.h5")
    example_itp = hdf5_to_interpolator(example_flow_field_file)
    example_field = StokesFlowField(example_itp, 293.15, 4.27e-26, 1.76e-5)
    example_dists = (
        x = Normal(0.0, 1e-3), z = Dirac(-0.128),
        vx = Normal(0.0, 0.1), vz = Dirac(10.0),
        r = Normal(13.5e-9, 2e-9), ρ = Dirac(1050.0)
    )
    example_source = SamplingSource{SphericalParticle2D{Float64}}(example_dists)
    example_detectors = Tuple(CMInject.SectionDetector{Float64,:z}.(-0.003:0.001:0.003, true))
    example_experiment = CMInject.Experiment(;
        source=example_source,
        n_particles=1000,
        fields=(example_field,),
        detectors=example_detectors,
        solver=EulerHeun(), time_span=(0.0, 0.05), time_step=1e-5, ensemble_alg=EnsembleThreads()
    )

    return (
        field = example_field,
        field_interpolator = example_itp,
        source = example_source,
        source_distributions = example_dists,
        detectors = example_detectors,
        experiment = example_experiment
    )
end

end
