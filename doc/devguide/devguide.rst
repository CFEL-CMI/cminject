.. _developer-guide:

CMInject Developer Guide
========================

This developer's guide gives a short overview of the steps required to develop CMInject and the
potential starting points for adding changes. For an overview of CMInject itself, please refer
to :ref:`index`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Installation
------------

To develop for CMInject, you should first install it as described in :ref:`user-guide`, making
sure to use ``python3 setup.py develop`` instead of ``python3 setup.py install``. This ensures
you don't have to reinstall after every code change you make to see an effect when importing or
running your changed code.

.. warning::
  When making changes to code that is compiled, like Cython code, you still need to rerun the
  ``python3 setup.py develop`` to see an effect.

Tests
-----

Tests are located in the ``tests/`` subdirectory. They can be executed
via the command ``pytest``. You may need to install ``pytest`` first,
as it is not an explicit dependency of CMInject:

.. code-block:: bash

    pip install pytest


You should then be able to just run ``pytest`` from your console. Make
sure you are in the top-level directory of the CMInject code folder
when executing ``pytest``, so that pytest can find the test files.


Starting points
---------------

You can find the full documentation of all modules in :ref:`the APIdoc<apidoc>`.

New concrete object definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Almost all of CMInject is written in terms of `abstract base classes (ABCs)
<https://docs.python.org/3/library/abc.html>`_. When defining a new kind of physical object like
a Field, Device, or Particle, please derive a subclass of the corresponding ABC. For a list of
all base classes, see :ref:`base-classes`.

For every ABC, there is documentation of concrete subclasses that you can refer to, to get a
better idea of what a specific kind of class does and how it can be implemented.

New setup definitions
~~~~~~~~~~~~~~~~~~~~~
After you implement all required abstract methods defined by the ABC, you can instantiate your
new subclass and use it. The next typical thing to define is a new kind of experimental setup,
based on the Setup subclass (see :ref:`list-of-setups`) and a combination of existing or newly
implemented subclasses of Particle, Field, etc.

Below, we've included verbatim the definition of :class:`cminject.setups.example.ExampleSetup`,
which is fairly simple. It defines a setup with a very small set of exposed parameters (``-f``,
``--rho`` and ``-pos``), consisting of only one ``Source``, one ``Device``, and three instances of
``Detector``.

This code is heavily annotated and should serve as a starting point for new setup definitions.
You can define a new setup by copying over this code to a differently named class. Then you can,
for example,

   * expose more parameters or remove existing ones(defined in ``get_parser()``, used in
     ``construct_experiment()``)
   * add more devices (see :mod:`cminject.devices`)
   * add more sources (see :mod:`cminject.sources`)
   * add instances of new object definitions that you've defined yourself. This can, e.g., be a new
     kind of ``Device``. Consider the relevant :ref:`base classes<base-classes>` and their abstract
     interfaces, which you must implement in your subclasses.

.. include:: ../../lib/cminject/setups/example.py
   :literal:


Performance considerations
--------------------------

A few points regarding runtime performance within this framework are discussed below.

Micro-optimizations
~~~~~~~~~~~~~~~~~~~
What is typically frowned upon and called an unnecessary micro-optimization can make or break the
performance of your simulation in a numerical simulation framework. Whether a method to calculate
a force takes 1us or 2us becomes highly relevant to the overall runtime of the program when this
method is called 10^8 times during the course of a simulation.

To achieve decent performance, you can time primitive calls and compare your possible
implementations. An overview of various useful profiling tools to help you with this task
`is given here
<https://jakevdp.github.io/PythonDataScienceHandbook/01.07-timing-and-profiling.html>`_.

Default simulation time step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Make sure to define an appropriate default time step for your Setup subclasses. Different
physical problem warrant different amounts of precision. While choosing a time step that is too
large will create numerical errors in user's simulation results, choosing one that is too small
can really impact performance negatively.

Particle instance size
~~~~~~~~~~~~~~~~~~~~~~
Make sure to not create Particle subclasses whose instances will have a large memory footprint,
unless absolutely necessary. Since the program is parallelized, Particle instances must be
serialized and deserialized between the main process and worker processes every time the
simulation of each single Particle starts and ends. Doing so requires I/O, which blocks your
thread so useful work can not be done during this (de)serialization step.


Contact
-------

Please direct any developer-centric questions to the authors,
`Simon Welker <mailto:simon.welker@cfel.de>`_ and
`Muhamed Amin <mailto:muhamed.amin@cfel.de>`_.
