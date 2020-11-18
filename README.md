# CMInject: A framework for particle injection trajectory simulations

CMInject is a Python framework that can be used to define and run particle trajectory simulations
for sample injection. It can run 1D, 2D, and 3D simulations based on the ODE describing Newton's
equation of motion, given that all parts in a simulation setup allow this dimensionality.

It is meant to be used by different groups of people and match different needs, and so it is meant
to be both:

- a library that provides a generic way to define and run whole simulation setups, and
- a collection of standard setups for immediate use.

The minimum required Python version to use it is 3.6.

CMInject is free software, see the [licensing conditions](LICENSE.md) for details.


## Installation

This section contains a short quickstart-style installation guide.
 
For detailed installation instructions and more, please generate the documentation (see section
_Generating docs_ below), and there, read the _Installation_ section in the _User Guide_.

### Prerequisites

CMInject requires the following packages before installation:

- Cython >= 0.29.10
- numpy >= 1.16.0

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

If you plan on developing CMInject you should use
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


<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 100
End:
-->
