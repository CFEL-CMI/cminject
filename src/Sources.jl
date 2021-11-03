using Distributions

abstract type AbstractSource end
abstract type AbstractSamplingSource <: AbstractSource end

struct SamplingSource{PT,Dists<:NamedTuple} <: AbstractSamplingSource where PT<:AbstractParticle
    distributions::Dists
end
SamplingSource{PT}(dists) where PT = SamplingSource{PT,typeof(dists)}(dists)

# Note that PT should either be StarkParticle or StarkParticle2d
struct StarkSamplingSource{PT, T} <: AbstractSamplingSource where {PT<:AbstractParticle, T}
    distributions::Dict{Symbol, Samp} where
        {Samp<:Sampleable{F, S} where {F<:VariateForm, S<:ValueSupport}}
    # Contains the distributions for the quantum states of the molecules
    stateDistributions::Dict{Symbol, Samp} where
        {Samp<:Sampleable{F, S} where {F<:VariateForm, S<:Discrete}}
    # Contains the particle properties which are the same for all particles
    particleProperties::Dict{Symbol, T}
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
    # TODO: Guarantee that PT is a Stark particle
    particles = [PT{typeof(starkCurves[1])}(; starkCurve=starkCurves[i], samples[i]...) for i=1:n]
end

function getStarkCurves(source::StarkSamplingSource{PT}, n::I) where {PT, I<:Integer}
    stateSamples = getSamples(source.stateDistributions, n)
    uniqueParameters = unique(stateSamples)
    starkCurves = Dict([(quantumParams,
                         # TODO: Parametric file name
                         interpolateStarkCurve("../cmistark/OCS.molecule"; quantumParams...))
                         for quantumParams ∈ uniqueParameters])
    [starkCurves[sample] for sample ∈ stateSamples]
end

function getSamples(dists::Dict{Symbol, Samp}, n::I) where {Symbol, Samp, I<:Integer}
    samples = [NamedTuple((sym, rand(dists[sym])) for sym in keys(dists)) for i=1:n]
end
