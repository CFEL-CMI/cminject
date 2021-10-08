using Distributions

abstract type Source end

struct SamplingSource{PT} <: Source where PT<:AbstractParticle
    distributions::Dict{Symbol, Samp} where {Samp<:Sampleable{F, S} where {F<:VariateForm,S<:ValueSupport}}
end

function generate(source::SamplingSource{PT}, n::I) where {PT, I<:Integer}
    dists = source.distributions
    samples = [NamedTuple((sym, rand(dists[sym])) for sym in keys(dists)) for i=1:n]
    particles = [PT(; samples[i]...) for i=1:n]
    particles
end