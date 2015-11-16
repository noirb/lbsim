# lbsim
A small &amp; simple Lattice-Boltzmann Method fluid simulator supporting complex boundaries. Produced as part of a lab course in Computational Fluid Dynamics at TUM from the [Chair of Scientific Computing](http://www5.in.tum.de/wiki/index.php/Home).

    Usage:
    
    lbsim <path_to_domain_description_file>

See the inputFiles directory for some examples of how to specify the domain parameters. 

There is also a graphical editor which can be used to create domain descriptions available in the [lb_ed repository](https://github.com/noirb/lb_ed).

Results of the simulation will be stored in (legacy) .vtk files in the same directory as the input file (prefixed with the input file's name), but there is also the possibility to visualize results as they're computed using the [FuLBLINKy project](https://github.com/noirb/FuLBLINKy).

## Build Instructions

The project should build without issues on Linux (or on Windows with Cygwin) by just issuing the `make` command from the project directory. Use `make debug` to include debugging symbols.
