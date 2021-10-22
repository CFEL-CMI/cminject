using DifferentialEquations

abstract type AbstractDetector
end

#= struct Hits{HitT,VecT<:AbstractVector{HitT}}
    __vec::VecT
end

struct DetectorHits{HitsT,VecT<:AbstractVector{HitsT}}
    __vec::VecT
end

function Base.getindex(hits::Hits{HitT,VecT}, idx) where {HitT,VecT<:AbstractVector{HitT}}
    [getindex(hit, idx) for hit in hits.__vec]
end
function Base.getindex(detector_hits::DetectorHits{HitsT,VecT}, idx) where {HitsT,VecT<:AbstractVector{HitsT}}
    [getindex(hits, idx) for hits in detector_hits.__vec]
end =#


function calculate_hits(
    solution::EnsembleSolution{A,B,Vector{RODESolution{C,D,Vector{Eltype},E,F,G,H,I,J,K,L}}},
    detector::Det
) where {Det<:AbstractDetector,Eltype,A,B,C,D,E,F,G,H,I,J,K,L}
    hits = Eltype[]  # initialize empty vector
    for one_solution ∈ solution
        trajectory = one_solution.u
        add_hits!(hits, trajectory, detector)
    end
    hits
end


struct SectionDetector{T,Dim} <: AbstractDetector
    at::T
    once::Bool
end

@generated function add_hits!(
    hits::AbstractVector{E}, trajectory::AbstractVector{E}, detector::SectionDetector{T,Dim}
) where {Dim, Syms, T, N, VT, E<:LArray{T,N,VT,Syms}}
    dim_idx = findfirst((sym) -> sym == Dim, Syms)
    if isnothing(dim_idx)
        error("Dimension $(Dim) does not exist in trajectory!")
        return :()
    end
    quote
        # Iterate over pairs of neighboring points (with a step of 1 point)
        @inbounds for (p₁, p₂) ∈ zip(trajectory[begin:end-1], trajectory[begin+1:end])
            found = false
            ps₁, ps₂ = p₁[$dim_idx], p₂[$dim_idx]
            if isnan(ps₁) || isnan(ps₂) return end # nothing to do anymore from here on
            Δps₁, Δps₂ = detector.at - ps₁, detector.at - ps₂

            # TODO perhaps rewrite this using better (more readable/reusable/accurate) interpolation
            if Δps₁ == 0 # p₁ is inside the section: just push directly
                push!(hits, p₁)
                found = true
            elseif Δps₁ < 0 && Δps₂ >= 0 || Δps₁ > 0 && Δps₂ <= 0
                # crossed detector going p₁ -> p₂ or p₂ -> p₁: linearly interpolate
                # swap points if the first is after the detector and the second is before/on the detector
                if Δps₁ < 0 && Δps₂ >= 0
                    ps₁, ps₂ = ps₂, ps₁
                end
                w₁ = (ps₂ - detector.at) / (ps₂ - ps₁)
                w₂ = 1 - w₁
                phase_space_point = w₁*p₁ + w₂*p₂
                phase_space_point[$dim_idx] = detector.at # overwrite to avoid tiny numerical offsets from interpolation
                push!(hits, phase_space_point)
                found = true
            end

            if found && detector.once
                break
            end
        end
    end
end

