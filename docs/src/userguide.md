# [User Guide](@id userguide)

## Setup

Please follow these steps:

- Install Julia 1.6 on your machine, [see here for instructions](https://julialang.org/downloads/)
- Clone or download [CMInject](https://github.com/CFEL-CMI/cminject)
- Open a terminal and enter the `cminject/` folder
- Run `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Enter the `build/` subfolder
- Execute `./create_sysimage.sh` -- this may take a while, but it will make your simulations much faster!
- You should now be able to run `./cminject.sh` without errors occurring. Use `./cminject.sh -h` to get help with running the program.

!!! tip "HDF5 errors"
    If you got some error message about hdf5, please install HDF5 headers and then try again.
    Possible commands to install HDF5 headers are:

    - `sudo apt-get install libhdf5-dev` on Ubuntu
    - `brew install hdf5` on Mac OS X

!!! note "Platform support"
    CMInject.jl currently supports Linux and Mac OS X.

## Running CMInject

Let's look at an example command to run a simulation with CMInject. It's a 2D simulation of spherical particles through a single ALS flow field:

```bash
./cminject.sh \
    -f data/nitrogen_1.8mbar_extended_nan.h5 \
    --x "G[0,1e-3]" --z -0.1285 --vx "G[0,0.1]" --vz "G[10,1]" \
    --r "G[13.5e-9,2e-9]" --rho 19320 \
    -d -0.003 -0.002 -0.001 0.0 0.001 0.002 0.003 \
    -t 0 0.05 --dt 1e-5 \
    -n 10000 \
    -o example_results.hdf5 \
    --plot
```
`-f` determines the flow-field file to use - in this case, an example ALS flow field delivered with this repository in the `data/` subfolder. The options `--x`, `--z`, `--vx`, `--vz` determine the position and velocity distributions: For example, `--x "G[0,1e-3]"` means that the `x` position will be distributed like a Gaussian with µ=0 and σ=0.001. `--z -0.1285` means that the `z` position will be exactly -0.1285 for all particles. `--r` and `--rho` determine the distribution of the particle radius and material density.

The `-d` option will put detectors at each given `z` position, in this case from -0.003 to 0.003 in steps of 0.001. `-t 0 0.05 --dt 1e-5` means that we will simulate for the timespan [0s, 0.05s] with an integrator time step of 10^-5.

`-n 10000` tells CMInject to simulate 10,000 particles, `-o example_results.hdf5` makes it write results (detector hits and particle properties) to the HDF5 file `example_results.hdf5`, and using `--plot` will show a plot of 100 particle trajectories and detector hits (this is of course optional).

Another example of running a Stark simulation:

```bash
./cminject.sh \
    -n 10000 \
    -D 3 \
    --x G[0,0.00052] \
    --y G[0,0.00052] \
    --z 0 \
    --vx G[0,10] \
    --vy G[0,10] \
    --vz G[1860,20] \
    --m 1.412436264e-25 \
    --J D[0,0] \
    --M D[0,0] \
    -t 0 2e-3 \
    --dt 1e-6 \
    -e test/example_field.h5 \
    -s test/pyrrole-water.molecule \
    -o results.h5 \
    --Boundaries S[0,0,0.065,0.0015,0.002] \
                 S[0,0,0.30172,0.00075,0.002] \
                 S[0,0,0.53613,0.00075,0.002] \
                 K[0,0.511,-1,0.002] \
                 C[-0.00198571,0.00198571,-0.00398571,0.00198571,0.34572,0.49972,0.002] \
    -d 0.71247 \
    --Solver RK4
```

Help can be obtained with the `--help` option.

## Running CMInject on Maxwell

There's an example script, `cminject_maxwell.sh`, to run CMInject on Maxwell. You should copy it to your home directory on Maxwell, so you can change any options you like.

You can submit a CMInject simulation job to Maxwell like so:

```bash
sbatch cminject_maxwell.sh \
    -f data/nitrogen_1.8mbar_extended_nan.h5 \
    --x "G[0,1e-3]" --z -0.1285 --vx "G[0,0.1]" --vz "G[10,1]" \
    --r "G[13.5e-9,2e-9]" --rho 19320 \
    -d -0.003 -0.002 -0.001 0.0 0.001 0.002 0.003 \
    -t 0 0.05 --dt 1e-5 \
    -n 10000 \
    -o example_results.hdf5
```

Which is exactly the same command as above in "Running CMInject", except for two minor differences:

- It uses `sbatch cminject_maxwell.sh` rather than `./cminject.sh`, to submit the job to Maxwell
- It doesn't use the `--plot` option, which doesn't make sense on a Maxwell compute node

You can of course adjust any simulation parameters you want here.

!!! tip "Default partition and duration"
    By default, the jobs submitted with this script will run on the `cfel` partition for a maximum time of 1 hour. You may want to change this if the `cfel` partition is full or you have very long computations or many simulations to run with one job.

!!! tip "Job status emails"
    The provided script won't send you emails about the job success or failure -- if you want to receive these, need to add the `#SBATCH --mail-type ALL` and `#SBATCH --mail-user youremail` options (with the correct email address).

!!! warning "Manual execution"
    If you execute CMInject manually on a Maxwell node without submitting a batch script, make sure that you have loaded the `julia` module, otherwise running `./cminject.sh` will fail. [See here for instructions](https://confluence.desy.de/display/IS/2020/06/24/Julia+on+Maxwell).
