#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "LBDefinitions.h"

// used for initializing flags array
typedef enum { 
    VARY_NONE = 0x0,
    VARY_X = 0x1,
    VARY_Y = 0x2,
    VARY_Z = 0x4
    } vary_flags;

#ifdef _WIN32 // we build as C++ on Windows, so need this to use vary_flags in bitwise operations
inline vary_flags operator|(vary_flags a, vary_flags b)
{
    return static_cast<vary_flags>(static_cast<int>(a) | static_cast<int>(b));
}
inline vary_flags& operator|=(vary_flags& a, vary_flags b)
{
    return a = a | b;
}
inline vary_flags operator&(vary_flags a, vary_flags b)
{
    return static_cast<vary_flags>(static_cast<int>(a) & static_cast<int>(b));
}
inline vary_flags& operator&=(vary_flags& a, vary_flags b)
{
    return a = a & b;
}
#endif

/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *xlength,                       /* reads domain's x-dimension. Parameter name: "xlength" */
    int *ylength,                       /* reads domain's y-dimension. Parameter name: "ylength" */
    int *zlength,                       /* reads domain's z-dimension. Parameter name: "zlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[]                        /* argv[1] shall contain the path to the config file */
);


/* initialises the particle distribution functions and the flagfield */
int initialiseFields(double *collideField, double *streamField, flag_data *flagField, int xlength, int ylength, int zlength, char* cellDataFile);

// sets a region of the flagField to a specific value
void setFlags(flag_data *flagField, cell_flag flag, int xstart, int ystart, int zstart, double* cell_parameters,  int xlength, int ylength, int zlength, vary_flags varying);

// performs some tests to ensure the provided domain configuration is valid & ok to simulate
int validateFlags(flag_data *flagField, int xlength, int ylength, int zlength);

void logInputError(char* problem, int line_number, char* line_content);
#endif

