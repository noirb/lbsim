# lbsim
A small &amp; simple Lattice-Boltzmann Method fluid simulator supporting complex boundaries. 

    Usage:
    
    lbsim <path_to_domain_description_file>

See the inputFiles directory for some examples of how to specify the domain parameters.

There is also a graphical editor which can be used to create domain descriptions available in the [lb_ed repository](https://github.com/noirb/lb_ed).

Results of the simulation will be stored in (legacy) .vtk files in the same directory as the input file (prefixed with the input file's name), but there is also the possibility to visualize results as they're computed using the [FuLBLINKy project](https://github.com/noirb/FuLBLINKy).

Produced with the CFD Lab at TUM from the [Chair of Scientific Computing](http://www5.in.tum.de/wiki/index.php/Home).

## Build Instructions

The project should build without issues on Linux (or on Windows with Cygwin) by just issuing the `make` command from the project directory. Use `make debug` to include debugging symbols.

If you want to build and run the project on Windows without Cygwin, see the [fulblinky](https://github.com/noirb/FuLBLINKy) project for a Visual Studuio 2015 project which will build it as either an EXE or DLL.
