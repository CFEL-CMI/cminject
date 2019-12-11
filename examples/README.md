# Example setups

This folder contains a simple example experimental setup, `simple_setup.SimpleSetup`, to run and to base own
Setup definitions on. Other similarly simple setups can be added to this folder; "real" and more complex setups should 
be put into lib/cminject/setups if they should be available to anyone using this library.

*Please read the README.md file in the root directory, and copy this directory to some place outside the
source directory before proceeding!*

### Running SimpleSetup
First, change directory to the copy of this directory you created, and (temporarily) add it to the PYTHONPATH,
to make it available to the main program, `cminject`:
    
    cd <copied_examples_directory>
    export PYTHONPATH="$PYTHONPATH:."
 
Now let's simply use the example 2D field provided to quickly run the setup with the main program:

    cminject -s simple_setup.SimpleSetup -f 2d_example_Field.h5 -n 100 -o example_output.h5 --loglevel info
    
If all went well, we can visualize the results:

    cminject_visualize -T example_output.h5 -H r,z
    
This command will show both a qualitative trajectory plot and a quantitative histogram plot of the detectors.
    
Let's further look at the available options and how to extend them.
    
### Getting help
To get more information about available options, just run `cminject -h`.
This will document all the available parameters for the main program and its default setup.
Let's now instead look at the help specific to `SimpleSetup`:

    $ cminject -s simple_setup.SimpleSetup -h
    
    [...]
    
    --------------------- Help for the setup specific parser: ----------------------
    optional arguments:
      -f str, --filename str
                            The filename of the flow field (HDF5).
      --rho float           The density of the particle material [kg/m^3].
      
There are currently only two parameters available, which is likely far too little for a flexible and interesting setup.
Adding further parameters can be done by adding arguments inside `SimpleSetup`'s `get_parser()` method. Refer to the 
[argparse documentation](https://docs.python.org/3.7/library/argparse.html).

The parsed values will then available inside `SimpleSetup`'s `construct_experiment()` method
on the `args` parameter by the name defined by the new parser argument, and can be further passed along to
Devices, PropertyUpdaters, etc.

## Converting flow fields to use
To use other flow fields, it is likely necessary to convert them first.
The `2d_example_field.h5` was created by converting `2d_example_field.txt` via:

    cminject_txt-to-hdf5 -i 2d_example_field.txt -o 2d_example_field.h5 -d 2 -fG N -ft 293.15 -m

This will convert the flow field to an HDF5 file, mirror it around the Z axis, and store information about the field
(`-fG N` defines that the gas in the fluid is nitrogen, and `-ft 293.15` defines its temperature to be 293.15K).