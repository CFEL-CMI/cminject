# CMInject: A framework for particle injection trajectory simulations

CMInject is a Python 3 framework that can be used to define and run particle trajectory simulations
for sample injection. It can run 1D, 2D and 3D simulations based on the ODE describing Newton's
equation of motion, given that all parts in a simulation setup allow this dimensionality.

It is meant to be used by different groups of people and match different needs, and so it is meant
to be both:

- a library that provides a generic way to define and run whole simulation setups, and
- a collection of such setups (coming soon).

The minimum required Python version to use it is 3.6.


## Installation

TODO: I suggest that this is given as simple `setup.py` commands here -- for people that see this on
github -- and that we have the `setup.py` and the `pip` versions in the installation guide. Whowver
decides to use anaconda can also figure out how to use that very specific infrastructure...


Below is a quickstart installation for Anaconda users, copied from CMInject's Sphinx documentation.
For other installation instructions and more, please generate the documentation (see section
_Generating docs_ below) and read the _Installation_ section in the _User Guide_.

### Prerequisites

CMInject requires the following packages: ...


### Installing CMInject:

Installation of CMInject proceeds through standard Python procedures, i.e., in the code directory
run
```
python setup.py install
```
A user (non-admin) installation can be preformed through
```
python setup.py install --user
```

If you plan on developing CMInject you can exploit
```
python setup.py develop --user
```


## Documentation

The documentation is available on [readthedocs](add link). It can also be generated from this code
through running
```
python setup.py build_sphinx
```
The generated HTML documentation can then be viewed by opening `./build/sphinx/html/index.html` in a
browser.


## Overview of an experiment Setup

To define a new experiment Setup, a user creates a subclass of `cminject.definitions.base.Setup`, and must then
implement the abstract methods defined by this class. In doing so, the user defines both a parser for experiment
parameters and a way to construct an `Experiment` instance from the parameters to this parser and common parameters to
the main program (defined in the program, `bin/cminject`). An `Experiment` instance is constructed with the following:

- A list of `Source`s (which generate `Particle`s)
- A list of `Device`s (which consist of several `Field`s and a `Boundary`)
- A list of `Detector`s (optional)
- A list of `PropertyUpdater`s (optional)
- Additional parameters like the `number_of_dimensions`, the `time_interval` and `time_step`, which are otherwise set
  to defaults which might not fit the setup's needs.

The user can then run the setup via `bin/cminject -s <setupcls>`, where `<setupcls>` is replaced with a fully qualified
Python import name for the Setup subclass they defined. All classes mentioned in the list above are described in more
detail in the following, and also have more extensive Python documentation (see section "Generating docs" above).


## Overview of abstract base definitions

These definitions are all abstract, i.e. the classes mentioned lack concrete implementations for the
tasks they're supposed to perform. There are also some concrete implementations for these classes
included in the framework; a user can either extend from there, or write their own implementations
from scratch.

Given below each abstract definition is a short list of example implementations included with the framework.

### Particle

A Particle is some physical entity that's able to move in space, and whose movement will be simulated
and possibly tracked.

It is an extensible data container and should be used to store all relevant simulated properties. On
a most basic level, the particles only store a phase space position (a 2n vector, where n is the
number of spatial dimensions of the simulated experiment) and a mass. This can easily be extended to
any number of properties by deriving a subclass.

> - `SphericalParticle`: A generic particle with a radius and density (rho) to calculate the mass with
> - `ThermallyConductiveSphericalParticle`: Derived from SphericalParticle, with additional thermal properties.

### Source

A Source is something capable of generating an initial distribution of Particles to run a simulation
with. There can be multiple Sources within an Experiment setup, but they must all generate the same
type (instances of the same subclass) of Particles.

> - `GaussianDistributedSphericalSource`: A source generating a gaussian distribution of
`SphericalParticle` instances or instances of its subclasses; x/y/z, vx/vy/vz and r are randomly
distributed according to a gaussian distribution with given μ and σ.

### Device
A Device is a combination of Fields and a Boundary, and is supposed to model a real-world
device in a real-world experimental setup.

### Boundary
A Boundary is some shape that encloses an area or volume (depending on experiment dimensionality).

It is mainly used to describe the shape of a Device, and can tell whether a Particle is inside itself
or not. When a Particle is within a Device's Z-boundary but is *not* considered inside of the Device's
Boundary, it is considered having hit the "wall" of the device and its simulation stops at that moment.

> - `SimpleZBoundary`: A Boundary enclosing a finite interval along the Z axis, but otherwise unbounded.
> - `CuboidBoundary`: A cuboid-shaped Boundary, i.e. defined in terms of a finite interval for every dimension.
> - `GridFieldBasedBoundary`: A Boundary that encompasses the cuboid that a grid interpolation field's coordinates enclose.

### Field
A Field is a description of an _acceleration field_, i.e. an object that can calculate the acceleration
of a Particle when it is inside the Field.

Fields are, apart from special PropertyUpdaters, the only objects that have an actual effect on Particles'
velocity and thus position. The acceleration "felt" by any Particle at any moment in time is the sum of
the accelerations of all Fields that have an effect at the Particle's position.

> - `StokesFluidFlowField`: An acceleration field based on an interpolated drag force from a grid-based fluid flow field.

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
> - `BrownianMotionPropertyUpater`: Adds Brownian motion to the particle flight path based on some fluid flow field's
macroscopic properties and the time step size of the integrator.

### ZBounded

A ZBounded object is something that is bounded in the Z direction (which is usually considered the last
dimension in the position vector, so it's X/Y/Z in 3D, and R/Z in 2D). Users don't need to use this
class directly, as it is just a "Mixin class". All of the following classes are subclasses of ZBounded:

> - `Device`
> - `Boundary`
> - `Field`
> - `Detector`




<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 100
End:
-->
