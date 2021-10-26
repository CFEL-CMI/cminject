using Distributions

abstract type Source end
abstract type AbstractSamplingSource <: Source end

struct SamplingSource{PT} <: AbstractSamplingSource where PT<:AbstractParticle
    distributions::Dict{Symbol, Samp} where {Samp<:Sampleable{F, S} where {F<:VariateForm,S<:ValueSupport}}
end

# Note that PT should either be StarkParticle or StarkParticle2d
struct StarkSamplingSource{PT} <: AbstractSamplingSource where PT<:AbstractParticle
    distributions::Dict{Symbol, Samp} where {Samp<:Sampleable{F, S} where {F<:VariateForm,S<:ValueSupport}}
end

function generate(source::SamplingSource{PT}, n::I) where {PT, I<:Integer}
    samples = getSamples(source, n)
    particles = [PT(; samples[i]...) for i=1:n]
    particles
end

function generate(source::StarkSamplingSource{PT}, n::I) where {PT, I<:Integer}
    samples = getSamples(source, n)
    # TODO: Calculate Stark curves using passed quantum numbers
    starkCurve = interpolateStarkCurve(50, -100, [-10, -5, 0, 5, 10])
    # TODO: Guarantee that PT is a Stark particle
    # TODO: This fails, because the type of ITP (for some reason) can't get implied
    particles = [PT(; starkCurve=starkCurve, samples[i]...) for i=1:n]
    particles
end

function getSamples(source::AbstractSamplingSource, n::I) where {I<:Integer}
    dists = source.distributions
    samples = [NamedTuple((sym, rand(dists[sym])) for sym in keys(dists)) for i=1:n]
end
