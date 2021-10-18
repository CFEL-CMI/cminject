using Parameters
using LinearAlgebra
using SplitApplyCombine

@with_kw struct StarkCurve{T<:Real}
    ΔE::T
    E_min::T
    # TODO: Think about interpolation
    energies::Vector{T}
end

"""
    calculateStarkCurves(ΔE, E_min, E_max, J_min, J_max, M, B, D, μ)

Calculates the Stark curves for a linear top molecule with completed shells.
This is done for multiple `J` quantum numbers at once (from `J_min` to `J_max`), but only for one `M` quantum number.
Note that while usually the rotational constants are given in Hertz, here Joules are used directly.

# Arguments
- `ΔE`: The field strength step size
- `E_min`: The starting field strength
- `E_max`: The ending field strength
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum quantum number
- `M::I`: The M quantum number
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `D`: The centrifugal distortion constant. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]
"""
function calculateStarkCurves(ΔE, E_min, E_max,
        # TODO: It might be cleaner to just pass the particle directly
        J_min::I, J_max::I, M::I, B, D, μ) where I <: Integer
    fieldJEnergy = [calculateEnergies(J_min, J_max, M, E, B, D, μ) for E ∈ E_min:ΔE:E_max]
    # We now have field -> J -> energy, but want J -> field -> energy
    jFieldEnergy = invert(fieldJEnergy)
    StarkCurve.(float(ΔE), float(E_min), jFieldEnergy)
end

"""
    calculateEnergies(J_min, , J_max, M, E, B, D, μ)

Calculates the energies for a lienar top molecule with completed shells.
This is done for the given `J` quantum number range, `M` quantum number and field strength.
The result is ordered by the `J` quantum number in ascending order (which is also the ascending order
of the energies themselves).
Note that while usually the rotational constants are given in Hertz, here Joules are used directly.

# Arguments
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum quantum number
- `M::I`: The M quantum number
- `E`: The field strength
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `D`: The centrifugal distortion constant. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]
"""
function calculateEnergies(J_min::I, J_max::I, M::I, E, B, D, μ) where I <: Integer
    hamiltonian = getHamiltonian(J_min, J_max, M, E, B, D, μ)
    sort(eigvals(hamiltonian))
end

"""
    getHamiltonian(J_min, J_max, M, E, B, D, μ)

Calculates the hamiltonian for a lienar top molecule with completed shells.
This is done for the given range of `J` quantum numbers, the `M` quantum number and the field strength.
Note that this is only one block of the whole block-diagonalized hamiltonian consisting of all the `M` quantum numbers.
Note that while usually the rotational constants are given in Hertz, here Joules are used directly.

# Arguments
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum quantum number
- `M::I`: The M quantum number
- `E`: The field strength
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `D`: The centrifugal distortion constant. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]
"""
function getHamiltonian(J_min::I, J_max::I, M::I, E, B, D, μ) where I <: Integer
    fieldFree = [B * J*(J+1) - D * (J*(J+1))^2 for J ∈ J_min:J_max]
    offDiagonal = [-μ * E * √((J+1)^2-M^2) / √((2*J+1) * (2*J+3)) for J ∈ J_min:J_max-1]
    fieldFree, offDiagonal = promote(fieldFree, offDiagonal)
    SymTridiagonal(fieldFree, offDiagonal)
end
