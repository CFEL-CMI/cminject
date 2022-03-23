using DifferentialEquations

"""
Abstract supertype for all particle detectors in CMInject.

Subtypes `MyDetector <: AbstractDetector` should implement:

- `add_hits!(hits::AbstractVector{E}, trajectory::AbstractVector{E}), detector::MyDetector)`

*Author:* Simon Welker
"""
abstract type AbstractDetector
end

"""
    add_hits!(hits, trajectory, detector)

Appends all hits that a single `trajectory` generates on a `detector` to `hits`. Modifies `hits` in-place, and
returns nothing. If no hits are detected, this function should be a no-op.

*Author:* Simon Welker
"""
function add_hits!(hits, trajectory, detector)
    error(
        "Missing implementation of add_hits!(hits::$(typeof(hits)), " *
        "trajectory::$(typeof(trajectory)), detector::$(typeof(detector)))"
    )
end

"""
    calculate_hits(solution::EnsembleSolution, detector)

Returns a vector of all hits on `detector` for all trajectories in a `solution::EnsembleSolution`.
Uses `add_hits!` internally to achieve this.

*Author:* Simon Welker
"""
function calculate_hits(
    solution::EnsembleSolution{A,B,Vector{C}},
    detector::Det
) where {Det<:AbstractDetector,A,B,C}
    if (size(solution)[1] == 0 || size(solution[1])[1] == 0)
        return []
    end
    Eltype = typeof(solution[1][1])
    hits = Eltype[]  # initialize empty vector
    for one_solution ∈ solution
        trajectory = one_solution.u
        add_hits!(hits, trajectory, detector)
    end
    hits
end


"""
    SectionDetector{T,Dim}(at, once)

A detector that intersects the phase-space with the hyperplane where `pos[Dim]`==`at` (where `pos` is a position in
this phase-space). `T` specifies the type of `Dim` in the phase-space, and `Dim` should typically be a Symbol (or Int).
If `once` is `true`, only the first hit of each particle is detected, and all subsequent hits of the same particle
are ignored.

# Examples

    `SectionDetector{Float64,:z}(0.01, true)`

The above line defines a detector at z=0.01, which detects all particles that cross the plane at that z position.
Because `true` was passed for the `once` parameter, only each particle's first hit will be registered.

*Author:* Simon Welker
"""
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

