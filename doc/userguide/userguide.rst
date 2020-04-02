.. _user-guide:

###################
CMInject User guide
###################

Welcome to the User guide for CMInject! CMInject is a Python3 framework for defining and running
nanoparticle trajectory simulations.

.. contents::

************
Installation
************

We recommend using `Anaconda`_, if possible, because of its easy setup and performance-optimized package distributions.

For Anaconda users
------------------
- Set up a new Anaconda environment for CMInject::

    conda create -n cminject python=3.8

- Activate the Anaconda environment::

    conda activate cminject

- Install Cython and numpy (required for installation)::

    conda install Cython numpy

- Install CMInject:

  - If you plan on developing CMInject (see `setup.py develop`_)::

      python setup.py develop

  - If you plan on only using CMInject::

      python setup.py install

For Python users without Anaconda
---------------------------------

- Create a virtual environment with `venv`_ in a directory of your choice::

    python3 -m venv /some/dir

- Activate the virtual environment (run the applicable **Command**):

  +------------+-----------------+------------------------------------+
  | Platform   | Shell           | Command                            |
  +============+=================+====================================+
  | POSIX      | bash/zsh        | source /some/dir/bin/activate      |
  |            +-----------------+------------------------------------+
  |            | fish            + . /some/dir/bin/activate.fish      |
  |            +-----------------+------------------------------------+
  |            | csh/tcsh        + source /some/dir/bin/activate.csh  |
  |            +-----------------+------------------------------------+
  |            | PowerShell Core + source /some/dir/bin/Activate.ps1  |
  +------------+-----------------+------------------------------------+
  | Windows    | cmd.exe         | some\\dir\\Scripts\\activate.bat   |
  |            +-----------------+------------------------------------+
  |            | PowerShell      | some\\dir\\Scripts\\Activate.ps1   |
  +------------+-----------------+------------------------------------+

- Install Cython and numpy (required for installation)::

    pip install Cython numpy

- Install CMInject:

  - If you plan on developing CMInject (see `setup.py develop`_)::

      python setup.py develop

  - If you plan on only using CMInject::

      python setup.py install

****************
Running CMInject
****************
``cminject`` is the program for running simulations. Let's look at an example call of it,
split into lines for readability::

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

You can also run different experiment setups. The call above is for the default of
``-s 'cminject.definitions.setups.OneFlowFieldSetup'``, which simulates particles moving through
exactly one flow field. This and other available setups are listed in :ref:`list-of-setups`.

.. note::
  Different setups can have different sets of parameters. To look at the parameters for a different
  setup, you can run, for example,
  ``cminject -s cminject.definitions.setups.desyatnikov_photophoresis -h``.

Other options exist and can be listed by running ``cminject -h``. The output file ``output.h5`` can
be viewed with ``cminject_visualize`` or further analyzed with ``cminject_analyze-asymmetry``, and
more virtual detectors can be inserted into the results file after simulation with
``cminject_reconstruct-detectors``.

.. note::
  If you want more information about how particles progress through your simulation, you can add the
  option ``--loglevel info``, or for even more verbose output, ``--loglevel debug``.

.. warning::
  ``cminject`` only accepts an HDF5 file as a flow field (i.e., the ``-f`` argument).
  See `cminject_txt-to-hdf5` for information on how to convert TXT files to such HDF5 files.

List of utility programs
------------------------
There are other programs to further process, analyze and visualize simulation results stored
by ``cminject``. This section gives a list of all these programs contained in CMInject and
describes each of them.

cminject_txt-to-hdf5
~~~~~~~~~~~~~~~~~~~~
``cminject_txt-to-hdf5`` was written to convert TXT files describing a field as a regular grid,
like flow field files, to HDF5 files. For example, the COMSOL Multiphysics software writes
out such TXT files. The reason this is useful is that large TXT files are very slow to read in in
comparison to HDF5 files.

To convert a file, run ``cminject_txt-to-hdf5 -i <infile.txt> -o <outfile.h5> -d <dimensions>``.
For convenience, you can store arbitrary attributes on the converted .h5 file that can be read
by CMInject's code, so you don't need to pass them when running the program. A typical set of such
attributes to store is ``-fG`` and ``-ft``, which store the gas type and temperature the field
was defined with.

.. warning::
  If the TXT file you are converting was generated for axisymmetric data, it might only contain
  entries for positive coordinates (e.g., the r in r/z coordinates). Since ``cminject`` does not
  know about this fact, particles might well cross into "negative r" and be considered 'lost'
  since they are, coordinate-wise, outside of the field. In this case, please use the ``-m`` option
  for ``cminject_txt-to-hdf5``, which mirrors the available data around the axis of symmetry and
  thus allows simulations to work as expected.

  Note that after doing this and running a simulation, you might want to work only with the
  absolute value of the simulated r positions, e.g.::

      r = resulthdf5['particles/0']['trajectory'][0]
      r = np.abs(r)

cminject_visualize
~~~~~~~~~~~~~~~~~~
``cminject_visualize`` visualizes result files. After you've run a simulation with
``cminject [...] -o resultfile.h5``, you can visualize this result file by running
``cminject_visualize``. There are currently two options for visualizing results available:

  - A trajectory visualization, which can be shown with ``-T`` and optionally configured through
    other parameters starting with ``-T``. It shows both trajectories as curves, and detectors
    as scatter plots::

        cminject_visualize
          resultfile.h5        # For resultfile.h5...
          -T                   # ...show trajectory plots...
          -Tn 30               # ...of 30 randomly sampled particles,
          -Tc                  # using color coding for velocities

    .. image:: img/vis2d_velcolor.png
    .. image:: img/vis3d.png
    .. image:: img/vis3d_velcolor.png

  - A detector histogram visualization (1D or 2D), which can be shown with ``-H x,y [x,y ...]``::

        # Show histograms for all stored detectors in resultfile.h5,
        # for a collection of dimension pairs to be shown as histograms together.
        # When one dimension has a constant value (e.g. z), a 1D histogram
        # will be shown, otherwise a 2D histogram will be shown.
        cminject_visualize resultfile.h5 -H x,y  x,z  y,z  x,vx  y,vy

    .. image:: img/vishist_r-z_r-vr.png

cminject_reconstruct-detectors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``cminject_reconstruct-detectors`` adds detectors at arbitrary z positions to an existing result
file. For this reconstruction to work, it's required that the given result file has stored the
trajectories; otherwise, there is nothing to reconstruct detectors from.

An example call is as follows::

    cminject_reconstruct-detectors
      resultfile.h5        # Reconstruct and add to resultfile.h5:
      --zs 0.01 0.0 -0.01  # At the z positions {0.01, 0, -0.01},
      --xis 1 2            # the properties stored in each trajectory
                           # at indices 1 and 2 (likely x and y),
      --zi 3               # assuming that z is stored at index 3.

.. note::
  The reconstructed detectors don't necessarily have the same shape as the detectors that were
  defined during the original simulation, so they are not stored with them, but instead under the
  key ``reconstructed_detectors``. Tools like ``cminject_visualize`` currently don't work with them,
  so analyses of the reconstructed data must be conducted manually.

cminject_analyze-asymmetry
~~~~~~~~~~~~~~~~~~~~~~~~~~
``cminject_analyze-asymmetry`` prints out information about the asymmetry of a 2D distribution at
each stored detector. The output format can either be nicely formatted text to be human-readable, or
CSV with the ``--csv`` parameter, for further data processing. An example call::

    cminject_analyze-asymmetry
       resultfile.h5   # Print the analysis results for resultfile.h5,
       --x 0 --y 1     # using the stored property at index 0 as the first
                       # dimension and the one at index 1 as the second.

which prints, for example, the following output::

    -------------------- Detector 0 --------------------
    α: 0.199
    e₀ = 6.473e-06	 e₁ = 9.693e-06
    θ₀ = -0.451π	 θ₁ = -0.951π
    μx = -1.658e-05	 μy = -3.031e-05

    -------------------- Detector 1 --------------------
    α: 0.934
    e₀ = 3.877e-07	 e₁ = 1.132e-05
    θ₀ = -0.523π	 θ₁ = 0.977π
    μx = -2.867e-05	 μy = -3.195e-04

This output can instead be printed as machine-readable CSV by passing the ``--csv`` flag parameter.

**********
References
**********
.. target-notes::

.. _`setup.py develop`: https://setuptools.readthedocs.io/en/latest/setuptools.html#develop-deploy-the-project-source-in-development-mode
.. _venv: https://docs.python.org/3/library/venv.html
.. _Anaconda: https://www.anaconda.com/distribution/