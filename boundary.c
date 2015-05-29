#include "boundary.h"
#include <stdio.h>
#include <math.h>

//macro for computing f_inv
#define FINV(x, y, z, i) (collideField[INDEXOF((x) + LATTICEVELOCITIES[(i)][0], (y) + LATTICEVELOCITIES[(i)][1], (z) + LATTICEVELOCITIES[(i)][2], 18 - (i))])

//macro for the calculation of the impact from the moving wall for treating moving wall boundary
#define WALL(density, i, wallVelocity) (2*LATTICEWEIGHTS[(i)]*(density)*(dot2(LATTICEVELOCITIES[(i)], wallVelocity) / ((C_S*C_S) )))

void treatBoundary(double *collideField, flag_data* flagField, int xlength, int ylength, int zlength){
  //loop over the boundaries

    int flag;
    double den; // the variable to store the density of the current cell, which is needed for treating moving wall boundary
    int ind;    // temporary variable for storing indexes

    const int maxIndex = NUMBER_OF_LATTICE_DIRECTIONS * (xlength+2) * (ylength+2) * (zlength+2);

    // Inner cells. Needed only in case there are some boundary cells inside the fluid. 
    // For the simple cavity case the inner conditionals are always false
    for (int i = 0; i <= xlength+1; i++)
    {
        for (int j = 0; j <= ylength+1; j++)
        {
            for(int k = 0; k <= zlength+1; k++)
            {
                flag = flagField[FINDEXOF(i, j, k)].flag;
                if (flag != FLUID)
                {
                    for (int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        ind = INDEXOF(i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2], 0);
                        if (ind < 0 || ind >= maxIndex)
                        {
                            continue;
                        }
                        else if (i + LATTICEVELOCITIES[l][0] > xlength + 1 ||
                                 j + LATTICEVELOCITIES[l][1] > ylength + 1 ||
                                 k + LATTICEVELOCITIES[l][2] > zlength + 1 ||
                                 i + LATTICEVELOCITIES[l][0] < 0 ||
                                 j + LATTICEVELOCITIES[l][1] < 0 ||
                                 k + LATTICEVELOCITIES[l][2] < 0 )
                        {
                            continue;
                        }
                        else if (flagField[FINDEXOF(i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2])].flag == FLUID)
                        {
                            if (flag == NO_SLIP)
                            {
                                collideField[INDEXOF(i, j, k, l)] = FINV(i, j, k, l);
                            }
                            else if (flag == MOVING_WALL)
                            {
                                ind = INDEXOF(i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2], 0);
                                computeDensity(collideField+ind, &den);
                                collideField[INDEXOF(i, j, k, l)] = FINV(i, j, k, l) + WALL(den, l, flagField[FINDEXOF(i, j, k)].parms);
                            }
                        }
                    }
                }
            }
        }
    }
}