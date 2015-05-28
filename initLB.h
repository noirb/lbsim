#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "LBDefinitions.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *xlength,                       /* reads domain's x-dimension. Parameter name: "xlength" */
    int *ylength,                       /* reads domain's y-dimension. Parameter name: "ylength" */
    int *zlength,                       /* reads domain's z-dimension. Parameter name: "zlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    char *cellDataPath,                 /* path to file containing cell descriptions. Parameter name: "cellDataPath" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[]                        /* argv[1] shall contain the path to the config file */
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(double *collideField, double *streamField,int *flagField, int xlength, int ylength, int zlength, char* cellDataFile);

#endif

