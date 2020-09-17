# Example setups

This folder contains a simple example experimental setup, `example_setup.ExampleSetup`, to run and to base own 
Setup definitions on. Other similarly simple setups can be added to this folder; "real" and more complex setups should  
be put into lib/cminject/setups, committed and submitted via a pull request if they should be available to everyone.

*Please copy this directory to some place outside the source directory before proceeding, so you won't make any
 changes in your cloned version of CMInject.*

### Running ExampleSetup
First, change directory to the copy of this directory you created, and (temporarily) add it to the PYTHONPATH,
to make its Python import path, "`example_setup.ExampleSetup`", available to `cminject`:
    
    cd <copied_examples_directory>
    export PYTHONPATH="$PYTHONPATH:."
 
Now let's just use the example 2D field provided to quickly run the setup with `cminject`, the main program:

    cminject -s simple_setup.ExampleSetup -f 2d_example_field.h5 -n 100 -o example_output.h5 -T  --loglevel info
    
We can then visualize the results:

    cminject_visualize -T example_output.h5 -H x,z
    
This command will show both a qualitative trajectory plot and a quantitative histogram plot of the detectors.
Let's further look at the available options and how to extend them.

### Getting help
To get more information about available options, just run `cminject -h`.
This will document all the available parameters for the main program and its default setup.
Let's now instead look at the help specific to `ExampleSetup`:

    $ cminject -s simple_setup.ExampleSetup -h
    
    [...]
    
    --------------------- Help for the setup specific parser: ----------------------
    optional arguments:
      -f str, --filename str
                            The filename of the flow field (HDF5).
      --rho float           The density of the particle material [kg/m^3].

There are currently only two parameters available, which is likely far too little for a flexible and interesting setup.
Adding further parameters can be done by adding arguments inside `ExampleSetup`'s `get_parser()` method. Refer to the 
[argparse documentation](https://docs.python.org/3.7/library/argparse.html).

The parsed values will then available inside `ExampleSetup`'s `construct_experiment()` method
on the `args` parameter by the name defined by the new parser argument, and can be further passed along to
Devices, Actions, etc.
