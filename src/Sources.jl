using Distributions

abstract type AbstractSource end
abstract type AbstractSamplingSource <: AbstractSource end

struct SamplingSource{PT,Dists<:NamedTuple} <: AbstractSamplingSource where PT<:AbstractParticle
    distributions::Dists
end
SamplingSource{PT}(dists) where PT = SamplingSource{PT,typeof(dists)}(dists)

struct StarkSamplingSource{PT, T} <: AbstractSamplingSource where {PT<:Union{StarkParticle,StarkParticle2D}, T}
    distributions::Dict{Symbol, Samp} where
        {Samp<:Sampleable{F, S} where {F<:VariateForm, S<:ValueSupport}}
    # Contains the distributions for the quantum states of the molecules
    stateDistributions::Dict{Symbol, Samp} where
        {Samp<:Sampleable{F, S} where {F<:VariateForm, S<:Discrete}}
    # The (relative) path to the filename of the molecule HDF5 file containg its Stark curves
    starkCurveFilename::AbstractString
end

function generate(source::SamplingSource{PT, Dists}, n::Integer) where {PT,Dists}
    samples = getSamples(source.distributions, n)
    particles = [PT(; samples[i]...) for i=1:n]
end

function generate(source::StarkSamplingSource{PT}, n::I)::Vector{PT} where {PT, I<:Integer}
    if n == 0
        return []
    end
    samples = getSamples(source.distributions, n)
    starkCurves = getStarkCurves(source, n)
    particles = [PT(; starkCurve=starkCurves[i], samples[i]...) for i=1:n]
end

function getStarkCurves(source::StarkSamplingSource{PT}, n::I) where {PT, I<:Integer}
    stateSamples = getSamples(source.stateDistributions, n)
    uniqueParameters = unique(stateSamples)
    # TODO: Potentially we don't want to load the HDF5 file for every parameter set
    starkCurves = Dict([(quantumParams,
                         interpolateStarkCurve(source.starkCurveFilename; quantumParams...))
                         for quantumParams ∈ uniqueParameters])
    [starkCurves[sample] for sample ∈ stateSamples]
end

function getSamples(dists::Union{Dict{Symbol, Samp}, NamedTuple}, n::I) where {Symbol, Samp, I<:Integer}
    samples = [NamedTuple((sym, rand(dists[sym])) for sym in keys(dists)) for i=1:n]
end
