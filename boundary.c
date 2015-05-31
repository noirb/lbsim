#include "boundary.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//macro for computing f_inv
#define FINV(x, y, z, i) (collideField[INDEXOF((x) + LATTICEVELOCITIES[(i)][0], (y) + LATTICEVELOCITIES[(i)][1], (z) + LATTICEVELOCITIES[(i)][2], 18 - (i))])

//macro for mirroring for FREE_SLIP
#define FMIRROR(x, y, z, i, mirror) (collideField[INDEXOF((x) + LATTICEVELOCITIES[(i)][0], (y) + LATTICEVELOCITIES[(i)][1], (z) + LATTICEVELOCITIES[(i)][2], (mirror))])

//macro for the calculation of the impact from the moving wall for treating moving wall boundary
#define WALL(density, i, wallVelocity) (2*LATTICEWEIGHTS[(i)]*(density)*(dot2(LATTICEVELOCITIES[(i)], wallVelocity) / ((C_S*C_S) )))

void treatBoundary(double *collideField, flag_data* flagField, int xlength, int ylength, int zlength){
    double density_ref = 1.0;

    int flag;
    double den; // the variable to store the density of the current cell, which is needed for treating moving wall boundary
    int ind;    // temporary variable for storing indexes

    double feq_inlet[19];
    double feq_temp[19];
    double velocity_temp[3];

    // reflected lattice directions from specific faces
    int nMirrorDirs = 5;
    const int mirrorsrc_16[5] = {14, 15, 16, 17, 18};
    const int mirrormap_16[5] = { 0,  1,  2,  3,  4};
    const int mirrorsrc_02[5] = { 0,  1,  2,  3,  4};
    const int mirrormap_02[5] = {14, 15, 16, 17, 18};
    const int mirrorsrc_12[5] = { 4, 11, 12, 13, 18};
    const int mirrormap_12[5] = { 0,  5,  6,  7, 14};
    const int mirrorsrc_06[5] = { 0,  5,  6,  7, 14};
    const int mirrormap_06[5] = { 4, 11, 12, 13, 18};
    const int mirrorsrc_10[5] = { 3,  7, 10, 13, 17};
    const int mirrormap_10[5] = { 1,  5,  8, 11, 15};
    const int mirrorsrc_08[5] = { 1,  5,  8, 11, 15};
    const int mirrormap_08[5] = { 3,  7, 10, 13, 17};

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
                            switch (flag)
                            {
                                case NO_SLIP:
                                    collideField[INDEXOF(i, j, k, l)] = FINV(i, j, k, l);
                                    break;
                                case MOVING_WALL:
                                    ind = INDEXOF(i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2], 0);
                                    computeDensity(collideField+ind, &den);
                                    collideField[INDEXOF(i, j, k, l)] = FINV(i, j, k, l) + WALL(den, l, flagField[FINDEXOF(i, j, k)].parms);
                                    break;
                                case FREE_SLIP:
                                    if (abs(LATTICEVELOCITIES[l][0] + LATTICEVELOCITIES[l][1] + LATTICEVELOCITIES[l][2]) == 1) // only mirror across faces, never at an angle
                                    {
                                        const int* mirrors;
                                        const int* sources;
                                        switch (l)
                                        {
                                            case 16:
                                                sources = mirrorsrc_16;
                                                mirrors = mirrormap_16;
                                                break;
                                            case  2:
                                                sources = mirrorsrc_02;
                                                mirrors = mirrormap_02;
                                                break;
                                            case 12:
                                                sources = mirrorsrc_12;
                                                mirrors = mirrormap_12;
                                                break;
                                            case  6:
                                                sources = mirrorsrc_06;
                                                mirrors = mirrormap_06;
                                                break;
                                            case 10:
                                                sources = mirrorsrc_10;
                                                mirrors = mirrormap_10;
                                                break;
                                            case  8:
                                                sources = mirrorsrc_08;
                                                mirrors = mirrormap_08;
                                                break;
                                        }

                                        for (int mirror = 0; mirror < nMirrorDirs; mirror++)
                                        {
                                            collideField[INDEXOF(i, j, k, sources[mirror])] = FMIRROR(i, j, k, sources[mirror], mirrors[mirror]);
                                        }
                                    }
                                    break;
                                case INFLOW:
                                    computeFeq(&density_ref, flagField[FINDEXOF(i, j, k)].parms, feq_inlet);
                                    collideField[INDEXOF(i, j, k, l)] = feq_inlet[l];
                                    break;
                                case OUTFLOW:
                                    ind = INDEXOF(i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2], 0);
                                    computeDensity(collideField+ind, &den);
                                    computeVelocity(&(collideField[ind]), &den, velocity_temp);
                                    computeFeq(&density_ref, velocity_temp, feq_temp);
                                    collideField[INDEXOF(i, j, k, l)] = feq_temp[18 - l] + feq_temp[l] - FINV(i, j, k, l);
                                    break;
                                case PRESSURE_IN:
                                    // TODO
                                    break;
                            }
                        }
                    }
                }
            }
        }
    }
}