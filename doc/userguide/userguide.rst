CMInject's User guide
=====================

Welcome to the User guide for CMInject! CMInject is a Python3 framework for defining and running
nanoparticle trajectory simulations.

Installation
------------

When using a virtual environment (e.g. with Anaconda), CMInject and its collection of executables can be
installed by running ``python3 setup.py install``.

.. warning::
  For running this ``python3 setup.py install``, you need to already have the ``Cython`` and
  ``numpy`` Python packages installed. Otherwise, use ``pip`` or ``conda`` to install them first.

.. note::
  If you are not in a virtual environment, and don't have write access to the global Python
  installation, you should use ``python3 setup.py install --user``.

.. note::
  After installation, you should have all executables (e.g. ``cminject``) available in your
  Terminal. If installation succeeded, but you still can't run ``cminject``, you should first try
  opening a new terminal, and if that doesn't help, you should investigate where the executables
  were installed (the installation process typically prints this information) and add this directory
  to your PATH environment variable.

.. note::
  If you plan on also developing for CMInject, it's best to run ``python setup.py develop``.
  Doing this will prevent you from having to run ``python setup.py install`` again after every
  change you make.

Executables overview
--------------------

The most important part of CMInject for you as a user is the collection of executables you can use to run simulations,
and process, analyze and visualize those simulations' results. This section gives a list of all executables contained in
CMInject and describes each of them.

cminject
~~~~~~~~
The ``cminject`` executable is the only tool for running simulations. Let's look at a minimal call of it, split into lines for readability::

    cminject -n 100          # Run the simulation for 100 particles.

      -D 3                   # Use a spatial dimensionality of 3.

      -f flowfield.h5        # Use the flow field file flowfield.h5. Must be a 3D flow field
                             # since we specified -D 3.

      -rho 1050 -r 50e-9     # Simulate particles with a density of 1050kg/m^3
                             # and a radius of 50nm

      -p "G 0.0 1e-3" 0 0    # Randomly generate initial particle positions, the first dimension
                             # (x) being normally (gaussian) distributed with mu = 0m and
                             # sigma = 1mm, and the others (y, z) being fixed at 0m.

      -v "G 0.0 1.0" 0 -10.0 # Randomly generate initial particle velocities, the first dimension
                             # being normally (gaussian) distributed with mu = 0m/s and
                             # sigma = 1m/s, the second (y) fixed at 0m/s, and the third (z)
                             # fixed at -10.0m/s.

      -d 0 -0.1              # Insert virtual detectors at 0m and -1cm

      -T                     # Track and store trajectories

      -B                     # Enable Brownian motion

      -o output.h5           # Write results to output.h5

Other options exist and can be listed by running ``cminject -h``. The output file ``output.h5`` can
be viewed with ``cminject_visualize`` or further analyzed with ``cminject_analyze-asymmetry``, and
more virtual detectors can be inserted into the results file after simulation with
``cminject_reconstruct-detectors``.

.. note::
  If you want more information about how particles progress through your simulation, you can add the
  option ``--loglevel info``, or for even more verbose output, ``--loglevel debug``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
