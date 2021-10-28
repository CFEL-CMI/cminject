using Distributions

abstract type Source end
abstract type AbstractSamplingSource <: Source end

struct SamplingSource{PT} <: AbstractSamplingSource where PT<:AbstractParticle
    distributions::Dict{Symbol, Samp} where {Samp<:Sampleable{F, S} where {F<:VariateForm,S<:ValueSupport}}
end

# Note that PT should either be StarkParticle or StarkParticle2d
struct StarkSamplingSource{PT, T} <: AbstractSamplingSource where {PT<:AbstractParticle, T}
    distributions::Dict{Symbol, Samp} where
        {Samp<:Sampleable{F, S} where {F<:VariateForm, S<:ValueSupport}}
    # Contains the distributions for the quantum states of the molecules
    stateDistributions::Dict{Symbol, Samp} where
        {Samp<:Sampleable{F, S} where {F<:VariateForm, S<:Discrete}}
    # Contains the particle properties which are the same for all particles
    particleProperties::Dict{Symbol, T}
    # Defines the symmetry of the to-be-created particles
    symmetry::Symmetry
end

function generate(source::SamplingSource{PT}, n::I) where {PT, I<:Integer}
    samples = getSamples(source.distributions, n)
    particles = [PT(; samples[i]...) for i=1:n]
end

function generate(source::StarkSamplingSource{PT}, n::I) where {PT, I<:Integer}
    if n == 0
        return []
    end
    samples = getSamples(source.distributions, n)
    # TODO: Make ΔE etc. changeable
    starkCurves = getStarkCurves(source, n, 0.1, 0, 10)
    # TODO: Guarantee that PT is a Stark particle
    particles = [PT{typeof(starkCurves[1])}(; starkCurve=starkCurves[i], samples[i]...) for i=1:n]
end

function getStarkCurves(source::StarkSamplingSource{PT}, n::I,
        ΔE, E_min, E_max) where {PT, I<:Integer}
    stateSamples = getSamples(source.stateDistributions, n)
    # The calculation of the Stark curves is done for multiple Js at once
    Js = [sample.J for sample ∈ stateSamples]
    J_min = minimum(Js)
    J_max = maximum(Js)
    withoutJ = [filter(p->first(p) != :J, pairs(sample)) for sample ∈ stateSamples]
    uniqueParameters = unique(withoutJ)
    # It's assumed here that for every possible combination of M, K, etc. all values of J
    #  are possible and hence the performance isn't decremented
    starkCurves = Dict([(quantumParams,
                         StarkEffect.calculateStarkCurves(ΔE, E_min, E_max,
                             source.symmetry,
                             namedTupleToDict((J_min=J_min, J_max=J_max, quantumParams...)),
                             source.particleProperties))
                        for quantumParams ∈ uniqueParameters])
    [starkCurves[withoutJ[i]][Js[i]-J_min+1] for i=1:n]
end

function namedTupleToDict(namedTuple::NamedTuple)::Dict
    Dict(zip(keys(namedTuple), namedTuple))
end

function getSamples(dists::Dict{Symbol, Samp}, n::I) where {Symbol, Samp, I<:Integer}
    samples = [NamedTuple((sym, rand(dists[sym])) for sym in keys(dists)) for i=1:n]
end
