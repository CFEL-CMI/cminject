using Parameters

@with_kw struct StarkCurve
    ΔE::Float64
    E_min::Float64
    energies::Vector{Float64}
end
