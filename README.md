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
For detailed installation instructions and more, please have a look at the
[User Guide of the official documentation](
https://cminject.readthedocs.io/en/latest/userguide/userguide.html#installation).


### Prerequisites

On some systems, in particular for Python version 3.6, it may be necessary to install the [HDF5
libraries](https://www.hdfgroup.org/downloads/hdf5/) before attempting to install CMInject.


### Installing CMInject:

Installation of CMInject proceeds through standard Python procedures, i.e., in the code directory
run
```
python setup.py install
```

A user (non-admin) installation can be performed through
```
python setup.py install --user
```

If you plan on developing CMInject you should use
```
python setup.py develop [--user]
```

## Documentation

The documentation is available [on readthedocs](https://cminject.readthedocs.io).
It can also be generated from this repository by running
```
python setup.py build_sphinx
```
The generated HTML documentation can then be viewed by opening `./build/sphinx/html/index.html` in a
browser.


## Tests

Tests are located in the `tests/` subdirectory. They can be executed via the command `pytest`. You may
need to install `pytest` first, as it is not an explicit dependency of CMInject:

```bash
pip install pytest
```

You should then be able to just run `pytest` from your console. Make sure you are in the top-level
directory of the CMInject code folder when executing `pytest`, so that pytest can find the test
files.


<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 100
End:
-->
