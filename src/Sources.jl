using Distributions

abstract type AbstractSource end

struct SamplingSource{PT,Dists<:NamedTuple} <: AbstractSource where PT<:AbstractParticle
    distributions::Dists
end
SamplingSource{PT}(dists) where PT = SamplingSource{PT,typeof(dists)}(dists)

function generate(source::SamplingSource{PT, Dists}, n::Integer) where {PT,Dists}
    dists = source.distributions
    samples = [NamedTuple((sym, rand(dists[sym])) for sym in keys(dists)) for i=1:n]
    particles = [PT(; samples[i]...) for i=1:n]
    particles
end