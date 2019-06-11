# CMInject: A framework for particle injection trajectory simulations

CMInject is a Python3 framework that can be used to define and run particle trajectory simulations
for sample injection. It runs both 2D and 3D simulations, and is in principle agnostic to the number of spatial
dimensions: It's up to the concrete implementations of all parts of the simulated experimental setup to be able
to deal with any number of spatial dimensions.

It is meant to be used by different groups of people and match different needs, and so it is meant to be both

- a library that provides a generic way to define and run whole simulation setups, and
- a collection of such setups (coming soon).


## Overview of an Experiment setup
A user creates an `Experiment` instance and passes the constructor lists of (i.e. one list of each)

- `Source`s (which generate `Particle`s)
- `Device`s (which consist of a `Field` and a `Boundary`)
- `Detector`s
- `PropertyUpdater`s

which are all described in more detail below. They can then call the Experiment's `.run()` method
to actually run the simulation and get back a list of all Particle objects that were initially generated,
modified as they were during the course of the simulated experiment.

## Overview of abstract base definitions

These definitions are all abstract, i.e. the classes mentioned lack concrete implementations
for the tasks they're supposed to perform. There are also some concrete implementations for these classes
included in the framework; a user can either extend from there, or write their own implementations
from scratch.

Given below each abstract definition is a short list of example implementations included with the framework.

### Particle
A Particle is some physical entity that's able to move in space, and whose movement will be simulated
and possibly tracked.

It is an extensible data container and should be used to store all relevant simulated properties. On a
most basic level, the particles only store a phase space position (a 2n vector, where n is the number of
spatial dimensions of the simulated experiment) and a mass. This can easily be extended to any number
of properties by deriving a subclass.

> - `SphericalParticle`: A generic particle with a radius and density (rho) to calculate the mass with
> - `ThermallyConductiveSphericalParticle`: Derived from SphericalParticle, with additional thermal properties.

### Source
A Source is something capable of generating an initial distribution of Particles to run a simulation with.
There can be multiple Sources within an Experiment setup, but they must all generate the same type
(instances of the same subclass) of Particles.

> - `GaussianSphericalSource`: A source generating a gaussian distribution of `SphericalParticle` instances or instances
of its subclasses; x/y/z, vx/vy/vz and r are randomly distributed according to a gaussian distribution with given μ and σ.

### Device
A Device is a combination of a Field and a Boundary, and is supposed to model a real-world
device in a real-world experimental setup.

### Boundary
A Boundary is some shape that encloses an area or volume (depending on experiment dimensionality).

It is mainly used to describe the shape of a Device, and can tell whether a Particle is inside itself
or not. When a Particle is within a Device's Z-boundary but is *not* considered inside of the Device's
Boundary, it is considered having hit the "wall" of the device and its simulation stops at that moment.

> - `SimpleZBoundary`: A Boundary enclosing a finite interval along the Z axis, but otherwise unbounded.
> - `CuboidBoundary`: A cuboid-shaped Boundary, i.e. defined in terms of a finite interval for every dimension.

### Field
A Field is a description of an _acceleration field_, i.e. an object that can calculate the acceleration
of a Particle when it is inside the Field.

Fields are the only objects that have an actual effect on Particles' velocity and thus position.
The acceleration "felt" by any Particle at any moment in time is the sum of the accelerations of all
Fields that have an effect at the Particle's position.

> - `FluidFlowField`: An acceleration field based on an interpolated drag force from a grid-based fluid flow field.

### Detector
A Detector is something that will try detecting any Particle hitting it and can calculate (interpolate)
the actual position where the hit occurred.

The full collection of the Particle's phase space position and properties will be stored when such a
Particle-Detector hit happens.

> - `SimpleZDetector`: A detector located at some Z position, detecting along the full X/Y plane.

### PropertyUpdater
A PropertyUpdater is something that can update any properties of a Particle instance handed to the
PropertyUpdater's `update` method.

This can involve arbitrarily complex computations based on the time and the particle's
position and properties, or also not involve any calculation at all
(e.g. `TrajectoryPropertyUpdater` just takes the time and the Particle's phase space position
and appends them to a trajectory list).

> - `TrajectoryPropertyUpdater`: Tracks the particle phase position along its whole flight path.

### ZBounded
A ZBounded object is something that is bounded in the Z direction (which is usually considered the last
dimension in the position vector, so it's X/Y/Z in 3D, and R/Z in 2D). Users don't need to use this
class directly, as it is just a "Mixin class". All of the following classes are subclasses of ZBounded:

> - `Device`
> - `Boundary`
> - `Field`
> - `Detector`