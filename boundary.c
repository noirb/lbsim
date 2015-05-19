#include "boundary.h"
//#include "LBDefinitions.h"
//#include "computeCellValues.h"
#include <stdio.h>
#include <math.h>

#define FINV(xlength, x, y, z, i) (collideField[INDEXOF((xlength), (x) + LATTICEVELOCITIES[(i)][0], (y) + LATTICEVELOCITIES[(i)][1], (z) + LATTICEVELOCITIES[(i)][2], 18 - (i))])

#define _WALL(density, i, wallVelocity) (2*LATTICEWEIGHTS[(i)]*(density)*(dot2(LATTICEVELOCITIES[(i)], wallVelocity) / ((C_S*C_S) * sqrt(doti(LATTICEVELOCITIES[(i)], LATTICEVELOCITIES[(i)]) ))) )
#define WALL(density, i, wallVelocity) _WALL(density, i, wallVelocity) 
//__extension__({wallDetails(density, i, wallVelocity); _WALL(density, i, wallVelocity); })
/*
void wallDetails(const double density, int i, const double* wallVelocity)
{
    printf("  WALL\n--------\n");
    printf("i: %d\nCi: (%d, %d, %d)\nWi: %f\nDensity: %f\nCi . U: %f\nResult: %f\n------\n\n", i, LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2], LATTICEWEIGHTS[i], density, dot2(LATTICEVELOCITIES[(i)], wallVelocity), _WALL(density, i, wallVelocity));
}
*/
void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  //loop over the boundaries

    int flag;
    double den;
    int ind;

    const int maxIndex = NUMBER_OF_LATTICE_DIRECTIONS*(xlength+2)*(xlength+2)*(xlength+2);

    //inner cells
    for (int i = 0; i <= xlength+1; i++)
    {
        for (int j = 0; j <= xlength+1; j++)
        {
            for(int k = 0; k <= xlength+1; k++)
            {
                flag = flagField[FINDEXOF(xlength, i, j, k)];
                if (flag != FLUID)
                {
                    for (int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        ind = INDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2], 0);
                        if (ind < 0 || ind >= maxIndex)
                        {
                            continue;
                        }
                        else if (i + LATTICEVELOCITIES[l][0] > xlength+1 ||
                                 j + LATTICEVELOCITIES[l][1] > xlength+1 ||
                                 k + LATTICEVELOCITIES[l][2] > xlength+1)
                        {
                            continue;
                        }
                        else if (flagField[FINDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j +  LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2])] == FLUID)
                        {
                            if (flag == NO_SLIP)
                            {
                                collideField[INDEXOF(xlength,i, j, k, l)] = FINV(xlength, i, j, k, l);
                            }
                            else if (flag == MOVING_WALL)
                            {
                                computeDensity(collideField+ind, &den);
                                collideField[INDEXOF(xlength, i, j, k, l)] = FINV(xlength, i, j, k, l) + WALL(den, l, wallVelocity);
                            }
                        }
                    }
                }
                
                if (collideField[INDEXOF(xlength, i, j, k, 0)] < 0)
                {
                    printf("[%d, %d, %d] bad: %f\n", i, j, k, collideField[INDEXOF(xlength, i, j, k, 0)]);
                }
            }
        }
    }
}
