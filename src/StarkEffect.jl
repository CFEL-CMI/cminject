using Parameters
using LinearAlgebra
using SplitApplyCombine

@with_kw struct StarkCurve{T}
    ΔE::T
    E_min::T
    # TODO: Think about interpolation
    energies::Vector{T}
end

# For linear top (recognizable by the quantum numbers)
# with completed shells (hence integer quantum numbers)
function calculateStarkCurves(ΔE::T, E_min::T, E_max::T,
        # TODO: It might be cleaner to just pass the particle directly
        J_min::Int, J_max::Int, M::Int, B::T, D::T, μ::T)::Vector{StarkCurve{T}} where T
    fieldJEnergy = [calculateEnergies(J_min, J_max, M, E, B, D, μ) for E ∈ E_min:ΔE:E_max]
    # We now have field -> J -> energy, but want J -> field -> energy
    jFieldEnergy = invert(fieldJEnergy)
    StarkCurve.(ΔE, E_min, jFieldEnergy)
end

function calculateEnergies(J_min::Int, J_max::Int, M::Int, E::T,
        B::T, D::T, μ::T)::Vector{T} where T
    hamiltonian = getHamiltonian(J_min, J_max, M, E, B, D, μ)
    sort(eigvals(hamiltonian))
end

# J_min and J_max denote the range of the J quantum number,
# M is the specific value of the M quantum number,
# E is the field strength,
# B is a rotational constant and
# D is the centrifugal distortion constant
function getHamiltonian(J_min::Int, J_max::Int, M::Int, E::T,
        B::T, D::T, μ::T)::Matrix{T} where T
    fieldFree = [B * J*(J+1) - D * (J*(J+1))^2 for J ∈ J_min:J_max]
    offDiagonal = [-μ * E * √((J+1)^2-M^2) / √((2*J+1) * (2*J+3)) for J ∈ J_min:J_max-1]
    Tridiagonal(offDiagonal, fieldFree, offDiagonal)
end
