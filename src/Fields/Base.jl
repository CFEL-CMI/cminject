"""
Abstract supertype for all acceleration fields in CMInject.

Subtypes `MyField <: Field` should implement:

- `acceleration(particle, field::MyField, time)`, for all applicable particle types

They may also implement, if there is a noise term:

- `noise(particle, field::MyField, time)`, for all applicable particle types

Note that both methods must return `NamedTuple`s, with the names of the affected parts of the phase-space position
as the keys. See the documentation of `acceleration` and `noise` for more details.
"""
abstract type Field
end

"""
    acceleration(particle::AbstractParticle, field::Field, time)

Calculates the acceleration that a `field` exerts on a `particle` at a given `time`.

Methods for this function must always return `NamedTuple`s, with the names of the affected quantities of the
phase-space position as the keys, and the change (differential) in these quantities as the associated values.

!!! warning
    Acceleration fields should typically only affect velocities in the phase-space position, since the
    acceleration is the change in velocity. You may find exceptions to this rule for specific circumstances, but
    typically when you return a NamedTuple like `(x = ..., y = ...)`, you should think again if that's really what
    you meant here.

!!! note
    For a field that's generic and may work on different particle types, you should omit the type of the
    `particle` parameter. If you do so, you permit the field to work together with all particles that have (at least)
    the properties that your field will affect: a form of duck typing. In other cases, you may intentionally restrict
    the type of particle. In this case, if possible, try to specify a supertype of all particle types that the field is
    applicable for, to not be unnecessarily restrictive.

!!! warning
    Remember to document the actual physical meaning (and units) of the keys in the returned `NamedTuple` for your
    concrete implementation! Otherwise, users may end up combining fields and particles that mean different things by
    the same names, which will lead to errors that are hard to track down.

# Implementation example

    struct GravitationalField <: Field
        strength::Float64
    end

    function acceleration(particle, field::GravitationalField, time)
        # The comma is necessary if returning only one value!
        # Without it, the return type would not be a NamedTuple.
        return (vz = -strength,)
    end

Here, our new `GravitationalField` affects the `vz` part of the phase-space position, and thereby implements a constant
acceleration in the `z` direction.
"""
function acceleration(particle, field::Field, time)
    error(
        "Missing implementation of acceleration(particle::$(typeof(particle)), field::$(typeof(field)), time)"
    )
end

"""
    noise(particle::AbstractParticle, field::Field, time)

Calculates the diagonal noise spectral intensities that a `field` exerts on a `particle` at a given `time`.

Methods for this function must always return `NamedTuple`s, with the names of the affected quantities of the
phase-space position as the keys, and the spectral intensity of the noise for these quantities as the associated values.

Note that the noise from acceleration fields should typically only affect velocities in the phase-space position
(see the discussion in the docs for `acceleration`.)

# Implementation example

    struct DecreasingNoiseField <: Field
        strength::Float64
    end

    function noise(particle, field::DecreasingNoiseField, time)
        factor = field.strength * exp(-time)
        return (vx = factor * abs(particle.x), vy = factor * abs(particle.y))
    end

Here, our new `DecreasingNoiseField` adds noise to `vx` and `vy`, with a strength that is dependent on the `x` and `y`
positions, respectively, and exponentially decreases over time.
"""
noise(particle, field::Field, time) = ()  # no noise term for fields required by default: return an empty tuple
