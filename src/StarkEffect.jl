# This whole file is heavily inspired by the CMIstark project, see https://doi.org/10.1016/j.cpc.2013.09.001
# TODO: Validate results
using Parameters
using LinearAlgebra
using SplitApplyCombine

"""
    StarkCurve

Represents a Stark curve, that is the energies over a certain range of field strengths.

This is sampled (not continuous) with a step size of `ΔE`. The start is `E_min`.

The field strengths are given in Volts/Meter [V/m] and the energies in Joule [J].
"""
# TODO: If the result is in Joule, why is the multiplication with h skipped
#  for the matrix elements then?
@with_kw struct StarkCurve{T<:Real}
    ΔE::T
    E_min::T
    # TODO: Store interpolated
    energies::Vector{T}
end

"""
    calculateStarkCurves(ΔE, E_min, E_max, J_min, J_max, M, B, D, μ)

Calculates the Stark curves for a linear top molecule with completed shells.
This is done for multiple `J` quantum numbers at once (from `J_min` to `J_max`),
but only for one `M` quantum number.

Note that while usually the rotational constants are given in Hertz,
here Joules are used directly.

# Arguments
- `ΔE`: The field strength step size. Unit: Volts/Meter [V/m]
- `E_min`: The starting field strength. Unit: Volts/Meter [V/m]
- `E_max`: The ending field strength. Unit: Volts/Meter [V/m]
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum J quantum number
- `M::I`: The M quantum number
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `D`: The centrifugal distortion constant. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]
"""
function calculateStarkCurves(ΔE, E_min, E_max,
        # TODO: It might be cleaner to just pass the particle directly
        J_min::I, J_max::I, M::I, B, D, μ)::AbstractVector{StarkCurve} where I <: Integer
    fieldJEnergy = [calculateEnergies(J_min, J_max, M, E, B, D, μ) for E ∈ E_min:ΔE:E_max]
    # We now have field -> J -> energy, but want J -> field -> energy
    jFieldEnergy = invert(fieldJEnergy)
    StarkCurve.(float(ΔE), float(E_min), jFieldEnergy)
end

"""
    calculateStarkCurves(ΔE, E_min, E_max, J_min, J_max, M, K, B, AC, Δ_J, Δ_JK, Δ_K, μ)

Calculates the Stark curves for a symmetric top molecule with completed shells.
This is done for multiple `J` quantum numbers at once (from `J_min` to `J_max`),
but only for one `M` and one `K` quantum number.

Note that while usually the rotational constants are given in Hertz,
here Joules are used directly.

# Arguments
- `ΔE`: The field strength step size. Unit: Volts/Meter [V/m]
- `E_min`: The starting field strength. Unit: Volts/Meter [V/m]
- `E_max`: The ending field strength. Unit: Volts/Meter [V/m]
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum J quantum number
- `K::I`: The K quantum number
- `M::I`: The M quantum number
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `AC`: Either the first (for prolate molecules) or third (for oblate molecules)
        rotational constant, A/C=h/(8π^2⋅I_a/b). Unit: Joule [J]
- `Δ_J`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `Δ_JK`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `Δ_K`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]

# Returns
Returns a vector of the calculated Stark curves
"""
function calculateStarkCurves(ΔE, E_min, E_max,
        # TODO: Once the particles are implemented, they should probably be passed in directly
        J_min::I, J_max::I, M::I, K::I, B, AC,
        Δ_J, Δ_JK, Δ_K, μ)::AbstractVector{StarkCurve} where I <: Integer
    fieldJEnergy = [calculateEnergies(J_min, J_max, M, K, E, B, AC, Δ_J, Δ_JK, Δ_K, μ)
                    for E ∈ E_min:ΔE:E_max]
    # We now have field -> J -> energy, but want J -> field -> energy
    jFieldEnergy = invert(fieldJEnergy)
    StarkCurve.(float(ΔE), float(E_min), jFieldEnergy)
end

"""
    calculateEnergies(J_min, J_max, M, E, B, D, μ)

Calculates the energies for a linear top molecule with completed shells.
This is done for the given `J` quantum number range, `M` quantum number and field strength.
The result is ordered by the `J` quantum number in ascending order
(which is also the ascending order of the energies themselves).

Note that while usually the rotational constants are given in Hertz,
here Joules are used directly.

# Arguments
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum J quantum number
- `M::I`: The M quantum number
- `E`: The field strength. Unit: Volts/Meter [V/m]
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `D`: The centrifugal distortion constant. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]

# Returns
Returns a vector of the calculated energies. Unit: Joule [J]
"""
function calculateEnergies(J_min::I, J_max::I, M::I, E,
        B, D, μ)::AbstractVector where I <: Integer
    hamiltonian = getHamiltonian(J_min, J_max, M, E, B, D, μ)
    sort(eigvals(hamiltonian))
end


"""
    calculateEnergies(J_min, J_max, M, K, E, B, AC, Δ_J, Δ_JK, Δ_K, μ)

Calculates the energies for a symmetric top molecule with completed shells.
This is done for the given `J` quantum number range, `M` quantum number,
`K` quantum number and field strength.
The result is ordered by the `J` quantum number in ascending order
(which is also the ascending order of the energies themselves).

Note that while usually the rotational constants are given in Hertz,
here Joules are used directly.

# Arguments
- `J_min::I`: The minimum J quantum number
- `J_max::I`: The maximum J quantum number
- `M::I`: The M quantum number
- `K::I`: The K quantum number
- `E`: The field strength. Unit: Volts/Meter [V/m]
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `AC`: Either the first (for prolate molecules) or third (for oblate molecules)
        rotational constant, A/C=h/(8π^2⋅I_a/b). Unit: Joule [J]
- `Δ_J`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `Δ_JK`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `Δ_K`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]

# Returns
Returns a vector of the calculated energies. Unit: Joule [J]
"""
function calculateEnergies(J_min::I, J_max::I, M::I, K::I, E, B, AC,
        Δ_J, Δ_JK, Δ_K, μ)::AbstractVector where I <: Integer
    hamiltonian = getHamiltonian(J_min, J_max, M, K, E, B, AC, Δ_J, Δ_JK, Δ_K, μ)
    sort(eigvals(hamiltonian))
end

"""
    getHamiltonian(J_min, J_max, M, E, B, D, μ)

Calculates the Hamiltonian for a lienar top molecule with completed shells.
This is done for the given range of `J` quantum numbers,
the `M` quantum number and the field strength.
Note that this is only one block of the whole
block-diagonalized Hamiltonian consisting of all the `M` quantum numbers.

Note that while usually the rotational constants are given in Hertz,
here Joules are used directly.

# Arguments
- `J_min::I`: The minimum J quantum number. Must be bigger or equal to `M`.
- `J_max::I`: The maximum J quantum number
- `M::I`: The M quantum number
- `E`: The field strength. Unit: Volts/Meter [V/m]
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `D`: The centrifugal distortion constant. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]

# Returns
Returns an AbstractMatrix (specifically a symmetrical tridiagonal matrix)
representing the Hamiltonian.
"""
function getHamiltonian(J_min::I, J_max::I, M::I, E, B, D, μ)::AbstractMatrix where I <: Integer
    fieldFree = [B * J*(J+1) - D * (J*(J+1))^2 for J ∈ J_min:J_max]
    stark = [-μ * E * √((J+1)^2-M^2) / √((2*J+1) * (2*J+3)) for J ∈ J_min:J_max-1]
    fieldFree, stark = promote(fieldFree, stark)
    SymTridiagonal(fieldFree, stark)
end

"""
    getHamiltonian(J_min, J_max, M, K, E, B, AC, Δ_J, Δ_JK, Δ_K, μ)

Calculates the Hamiltonian for a symmetric top molecule with completed shells.
This is done for the given range of `J` quantum numbers,
the `M` quantum number, the `K` quantum number and the field strength.

Note that this is only one block of the whole
block-diagonalized Hamiltonian consisting of all the `M` and `K` quantum numbers.

Note that while usually the rotational constants are given in Hertz,
here Joules are used directly.

# Arguments
- `J_min::I`: The minimum J quantum number. Must be bigger or equal to `M` and `|K|`.
- `J_max::I`: The maximum J quantum number
- `M::I`: The M quantum number
- `K::I`: The K quantum number
- `E`: The field strength. Unit: Volts/Meter [V/m]
- `B`: The second rotational constant, B=h/(8π^2⋅I_b). Unit: Joule [J]
- `AC`: Either the first (for prolate molecules) or third (for oblate molecules)
        rotational constant, A/C=h/(8π^2⋅I_a/b). Unit: Joule [J]
- `Δ_J`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `Δ_JK`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `Δ_K`: One of the first order quartic centrifugal distortion constants. Unit: Joule [J]
- `μ`: The dipole moment along its symmetry axis. Unit: Coulomb⋅Meter [Cm]

# Returns
Returns an AbstractMatrix (specifically a symmetrical tridiagonal matrix) representing the Hamiltonian.
"""
function getHamiltonian(J_min::I, J_max::I, M::I, K::I, E, B, AC,
        Δ_J, Δ_JK, Δ_K, μ)::AbstractMatrix where I <: Integer
    rigid = [B * J*(J+1) + (AC-B) * K^2 for J ∈ J_min:J_max]
    distortion = [-Δ_J * (J*(J+1))^2 - Δ_JK * J*(J+1) * K^2 - Δ_K * K^4 for J ∈ J_min:J_max]
    fieldFree = rigid + distortion
    starkMain = [-M * K / (J*(J+1)) * μ * E for J ∈ J_min:J_max]
    starkOff = [-√((J+1)^2 - K^2) * √((J+1)^2 - M^2) / ((J+1) * √((2*J+1) * (2*J+3))) * μ * E
                for J ∈ J_min:J_max]
    fieldFree, starkMain, starkOff = promote(fieldFree, starkMain, starkOff)
    SymTridiagonal(fieldFree + starkMain, starkOff)
end
