using Parameters
using LinearAlgebra

@with_kw struct StarkCurve{T}
    ΔE::T
    E_min::T
    energies::Vector{T}
end

# For linear top (recognizable by the quantum numbers)
# with completed shells (hence integer quantum numbers)
function calculateStarkCurve(ΔE::T, E_min::T, E_max::T,
        J_min::Int, J_max::Int, M::Int)::Vector{StarkCurve{T}} where T
    # TODO
    [StarkCurve(ΔE=ΔE, E_min=E_min, energies=[1.0,2.0,3.0])]
end

function getHamiltonian(J_min::Int, J_max::Int, M::Int, E::T)::Matrix{T} where T
    # TODO
    Tridiagonal(collect(J_min+1:J_max), collect(J_min:J_max), collect(J_min:J_max-1))
end
