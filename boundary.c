#include "boundary.h"
#include <stdio.h>
#include <math.h>

//macro for computing f_inv
#define FINV(xlength, x, y, z, i) (collideField[INDEXOF((xlength), (x) + LATTICEVELOCITIES[(i)][0], (y) + LATTICEVELOCITIES[(i)][1], (z) + LATTICEVELOCITIES[(i)][2], 18 - (i))])

//macro for the calculation of the impact from the moving wall for treating moving wall boundary
#define _WALL(density, i, wallVelocity) (2*LATTICEWEIGHTS[(i)]*(density)*(dot2(LATTICEVELOCITIES[(i)], wallVelocity) / ((C_S*C_S) )))
#define WALL(density, i, wallVelocity) _WALL(density, i, wallVelocity) 

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  //loop over the boundaries

    int flag;
    double den; // the variable to store the density of the current cell, which is needed for treating moving wall boundary
    int ind; // temporary variable for storing indexes

    int N = 5;//number of fluid cells, which can be accessed from the domain boundary (sides of the cube)
    //using lattice directions
    const int xyo[5] = {14, 15, 16, 17, 18}; // Directions, allowed for the plane xy0 (z=0)
    const int xyn[5] = {0, 1 ,2, 3, 4}; // Directions, allowed for the plane xyn (z=n)
    const int xoz[5] = {4, 11, 12, 13, 18}; // Directions, allowed for the plane x0z (y = 0)
    const int xnz[5] = {0, 5, 6, 7, 14}; // Directions, allowed for the plane xnz (y=n)
    const int oyz[5] = {3, 7, 10, 13, 17}; // Directions, allowed for the plane xnz (z=0)
    const int nyz[5] = {1, 5, 8, 11, 15}; // Directions, allowed for the plane xnz (z=n)


    //inner cells. Needed only in case there are some boundary cells inside the fluid. In current assignment
    //both if-statements are always false
    for (int i = 1; i <= xlength; i++)
    {
        for (int j = 1; j <= xlength; j++)
        {
            for(int k = 1; k <= xlength; k++)
            {
                flag = flagField[FINDEXOF(xlength, i, j, k)];
                if (flag == NO_SLIP)
                {
                    for ( int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        collideField[INDEXOF(xlength,i, j, k, l)] = FINV(xlength, i, j, k, l);
                    }
                }
                else if (flag == MOVING_WALL)
                {
                    for (int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        ind = INDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j + LATTICEVELOCITIES[l][1], k + LATTICEVELOCITIES[l][2], 0);
                        computeDensity(collideField+ind, &den);
                        collideField[INDEXOF(xlength, i, j, k, l)] = FINV(xlength, i, j, k, l) + WALL(den, l, wallVelocity);
                    }
                }
            }
        }
    }

 // loop over different planes (sides of the cube) to update only the directions which point to the fluid
    //loop over xyo and XYN planes
    for(int i = 1; i <= xlength; i++)
    {
        for (int j = 1; j <= xlength; j++)
        {
            flag = flagField[FINDEXOF(xlength, i, j, 0)];
            if (flag == NO_SLIP) //update 14, 15, 16, 17, 18
            {
                for (int l = 0; l<N; l++)
                {
                    collideField[INDEXOF(xlength, i, j, 0, xyo[l])] = FINV(xlength, i, j, 0, xyo[l]);
                }
            }
            else if (flag == MOVING_WALL)
            {
                for(int l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xyo[l]][0], j + LATTICEVELOCITIES[xyo[l]][1], 0 + LATTICEVELOCITIES[xyo[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength, i, j, 0, xyo[l])] = FINV(xlength, i, j, 0, xyo[l]) + WALL(den, xyo[l], wallVelocity);
                }
            }

            flag = flagField[FINDEXOF(xlength, i, j, xlength+1)];
            if (flag == NO_SLIP) //update 0, 1, 2, 3, 4
            {
                for (int l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength, i, j, xlength+1, xyn[l])] = FINV(xlength, i, j, xlength+1, xyn[l]);
                }
            }
            else if (flag == MOVING_WALL)
            {
                for(int l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xyn[l]][0], j + LATTICEVELOCITIES[xyn[l]][1], xlength+1 + LATTICEVELOCITIES[xyn[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength, i, j, xlength+1, xyn[l])] = FINV(xlength, i, j, xlength+1, xyn[l]) + WALL(den, xyn[l], wallVelocity);
                }
            }
        }
    }

    //loop over OYZ and NYZ planes
    for(int j = 1; j <= xlength; j++)
    {
        for (int k = 1; k <= xlength; k++)
        {
            //0YZ
            flag = flagField[FINDEXOF(xlength, 0, j, k)];
            if (flag == NO_SLIP) //update 3, 7, 10, 13, 17
            {
                for (int l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength, 0, j, k, oyz[l])] = FINV(xlength, 0, j, k, oyz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 3, 7, 10, 13, 17
            {
                for(int l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, 0 + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength, 0, j, k, oyz[l])] = FINV(xlength, 0, j, k, oyz[l]) + WALL(den, oyz[l], wallVelocity);
                }
            }

            //NYZ
            flag = flagField[FINDEXOF(xlength, xlength+1, j, k)];
            if (flag == NO_SLIP) //update 1, 5, 8, 11, 15
            {
                for (int l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength, xlength+1, j, k, nyz[l])] = FINV(xlength, xlength+1, j, k, nyz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 1, 5, 8, 11, 15
            {
                for(int l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[nyz[l]][0], j+LATTICEVELOCITIES[nyz[l]][1], k+LATTICEVELOCITIES[nyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength, xlength+1, j, k, nyz[l])] = FINV(xlength, xlength+1, j, k, nyz[l]) + WALL(den, nyz[l], wallVelocity);
                }
            }
        }
    }

    //loop over X0Z and XNZ planes
    for(int i = 1; i <= xlength; i++)
    {
        for (int k = 1; k <= xlength; k++)
        {
            //X0Z
            flag = flagField[FINDEXOF(xlength, i, 0, k)];
            if (flag == NO_SLIP) //update 4, 11, 12, 13, 18
            {
                for (int l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength, i, 0, k, xoz[l])] = FINV(xlength, i, 0, k, xoz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 4, 11, 12, 13, 18
            {
                for(int l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i+LATTICEVELOCITIES[xoz[l]][0], 0+LATTICEVELOCITIES[xoz[l]][1], k+LATTICEVELOCITIES[xoz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength, i, 0, k, xoz[l])] = FINV(xlength, i, 0, k, xoz[l]) + WALL(den, xoz[l], wallVelocity);
                }
            }

            //XNZ
            flag = flagField[FINDEXOF(xlength, i, xlength+1, k)];
            if (flag == NO_SLIP) //update 0, 5, 7, 6, 14
            {
                for (int l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength, i, xlength+1, k, xnz[l])] = FINV(xlength, i, xlength+1, k, xnz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 0, 5, 7, 6, 14
            {
                for(int l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xnz[l]][0], xlength+1 + LATTICEVELOCITIES[xnz[l]][1], k+LATTICEVELOCITIES[xnz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength, i, xlength+1, k, xnz[l])] = FINV(xlength, i, xlength+1, k, xnz[l]) + WALL(den, xnz[l], wallVelocity);
                }
            }
        }
    }

//Loop over the edges of the cube, to update only one direction pointing to the fluid.
        /* ----- */
        /* Edges */
        /* ----- */

    //X-edges
    for(int i = 1; i <= xlength; i++)
    {
        //(i, 0, 0)
        flag = flagField[FINDEXOF(xlength, i, 0, 0)];
        if (flag == NO_SLIP) //update 18
        {
            collideField[INDEXOF(xlength, i, 0, 0, 18)] = FINV(xlength, i, 0, 0, 18);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[18][0], LATTICEVELOCITIES[18][1], LATTICEVELOCITIES[18][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, i, 0, 0, 18)] = FINV(xlength, i, 0, 0, 18) + WALL(den, 18, wallVelocity);
        }

        //(i, 0, xlength+1)
        flag = flagField[FINDEXOF(xlength, i, 0, xlength+1)];
        if (flag == NO_SLIP) //update 4
        {
            collideField[INDEXOF(xlength, i, 0, xlength+1, 4)] = FINV(xlength, i, 0, xlength+1, 4);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[4][0], LATTICEVELOCITIES[4][1], xlength+1 + LATTICEVELOCITIES[4][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, i, 0, xlength+1, 4)] = FINV(xlength, i, 0, xlength+1, 4) + WALL(den, 4, wallVelocity);
        }

        //(i, xlength+1, 0)
        flag = flagField[FINDEXOF(xlength, i, xlength+1, 0)];
        if (flag == NO_SLIP) //update 14
        {
            collideField[INDEXOF(xlength, i, xlength+1, 0, 14)] = FINV(xlength, i, xlength+1, 0, 14);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[14][0], xlength+1+LATTICEVELOCITIES[14][1], LATTICEVELOCITIES[14][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, i, xlength+1, 0, 14)] = FINV(xlength, i, xlength+1, 0, 14) + WALL(den, 14, wallVelocity);
        }

        //(i, xlength+1, xlength+1)
        flag = flagField[FINDEXOF(xlength, i, xlength+1, xlength+1)];
        if (flag == NO_SLIP) //update 0
        {
            collideField[INDEXOF(xlength, i, xlength+1, xlength+1, 0)] = FINV(xlength, i, xlength+1, xlength+1, 0);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[0][0], xlength+1+LATTICEVELOCITIES[0][1], xlength+1+LATTICEVELOCITIES[0][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, i, xlength+1, xlength+1, 0)] = FINV(xlength, i, xlength+1, xlength+1, 0) + WALL(den, 0, wallVelocity);
        }
    }


    //Y-edges
    for(int j = 1; j <= xlength; j++)
    {
        //(0, j, 0)
        flag = flagField[FINDEXOF(xlength, 0, j, 0)];
        if (flag == NO_SLIP) //update 17
        {
            collideField[INDEXOF(xlength, 0, j, 0, 17)] = FINV(xlength, 0, j, 0, 17);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[17][0], j+LATTICEVELOCITIES[17][1], LATTICEVELOCITIES[17][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, 0, j, 0, 17)] = FINV(xlength, 0, j, 0, 17) + WALL(den, 17, wallVelocity);
        }

        //(xlength+1, j, 0)
        flag = flagField[FINDEXOF(xlength, xlength+1, j, 0)];
        if (flag == NO_SLIP) //update 15
        {
            collideField[INDEXOF(xlength, xlength+1, j, 0, 15)] = FINV(xlength, xlength+1, j, 0, 15);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[15][0], j+LATTICEVELOCITIES[15][1], LATTICEVELOCITIES[15][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, xlength+1, j, 0, 15)] = FINV(xlength, xlength+1, j, 0, 15) + WALL(den, 15, wallVelocity);
        }

        //(0, j, xlength+1)
        flag = flagField[FINDEXOF(xlength, 0, j, xlength+1)];
        if (flag == NO_SLIP) //update 3
        {
            collideField[INDEXOF(xlength, 0, j, xlength+1, 3)] = FINV(xlength, 0, j, xlength+1, 3);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[3][0], j+LATTICEVELOCITIES[3][1], xlength+1+LATTICEVELOCITIES[3][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, 0, j, xlength+1, 3)] = FINV(xlength, 0, j, xlength+1, 3) + WALL(den, 3, wallVelocity);
        }

        //(xlength+1, j, xlength+1)
        flag = flagField[FINDEXOF(xlength, xlength+1, j, xlength+1)];
        if (flag == NO_SLIP) //update 1
        {
            collideField[INDEXOF(xlength,xlength+1, j, xlength+1, 1)] = FINV(xlength, xlength+1, j, xlength+1, 1);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[1][0], j+LATTICEVELOCITIES[1][1], xlength+1+LATTICEVELOCITIES[1][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, j, xlength+1, 1)] = FINV(xlength, xlength+1, j, xlength+1, 1) + WALL(den, 1, wallVelocity);
        }
    }


    //Z-edges
    for(int k = 1; k <= xlength; k++)
    {
        //(0, 0, k)
        flag = flagField[FINDEXOF(xlength, 0, 0, k)];
        if (flag == NO_SLIP) //update 13
        {
            collideField[INDEXOF(xlength, 0, 0, k, 13)] = FINV(xlength, 0, 0, k, 13);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[13][0], LATTICEVELOCITIES[13][1], k+LATTICEVELOCITIES[13][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength, 0, 0, k, 13)] = FINV(xlength, 0, 0, k, 13) + WALL(den, 13, wallVelocity);
        }

        //(xlength+1, 0, k)
        flag = flagField[FINDEXOF(xlength, xlength+1, 0, k)];
        if (flag == NO_SLIP) //update 11
        {
            collideField[INDEXOF(xlength,xlength+1, 0, k, 11)] = FINV(xlength, xlength+1, 0, k, 11);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[11][0], LATTICEVELOCITIES[11][1], k+LATTICEVELOCITIES[11][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, 0, k, 11)] = FINV(xlength, xlength+1, 0, k, 11) + WALL(den, 11, wallVelocity);
        }

        //(0, xlength+1, k)
        flag = flagField[FINDEXOF(xlength, 0, xlength+1, k)];
        if (flag == NO_SLIP) //update 7
        {
            collideField[INDEXOF(xlength,0, xlength+1, k, 7)] = FINV(xlength, 0, xlength+1, k, 7);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[7][0], xlength+1+LATTICEVELOCITIES[7][1], k+LATTICEVELOCITIES[7][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, xlength+1, k, 7)] = FINV(xlength, 0, xlength+1, k, 7) + WALL(den, 7, wallVelocity);
        }

        //(xlength+1, xlength+1, k)
        flag = flagField[FINDEXOF(xlength, xlength+1, xlength+1, k)];
        if (flag == NO_SLIP) //update 5
        {
            collideField[INDEXOF(xlength,xlength+1, xlength+1, k, 5)] = FINV(xlength, xlength+1, xlength+1, k, 5);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1+LATTICEVELOCITIES[5][0], xlength+1+LATTICEVELOCITIES[5][1], k+LATTICEVELOCITIES[5][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, xlength+1, k, 5)] = FINV(xlength, xlength+1, xlength+1, k, 5) + WALL(den, 5, wallVelocity);
        }
    }

}
