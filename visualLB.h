#ifndef _VISUALLB_H_
#define _VISUALLB_H_
#include <stdio.h>
#include "LBDefinitions.h"

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtkOutput(const double * const collideField, const flag_data * const flagField, const char * filename, unsigned int t, int xlength, int ylength, int zlength);

/**
 * Method for writing header information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkHeader( FILE *fp, int xlength, int ylength, int zlength);

/**
 * Method for writing grid coordinate information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkPointCoordinates( FILE *fp, int xlength, int ylength, int zlength, double dx, double dy, double dz);

#endif

