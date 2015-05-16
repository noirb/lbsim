#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

#define FINV(xlength, x, y, z, i) collideField[INDEXOF((xlength), (x) + LATTICEVELOCITIES[i][0], (y) + LATTICEVELOCITIES[i][1], (z) + LATTICEVELOCITIES[i][2], 18 - (i))]
#define WALL(density, i) 2*LATTICEWEIGHTS[(i)]*(density)*(dot2(LATTICEVELOCITIES[(i)], wallVelocity) / (C_S*C_S))

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  //loop over the boundaries

    int i,j, k, l;
    int flag;
    double den;
    int ind;
    int N = 5;
    const int xyo[5] = {14, 15, 16, 17, 18};
    const int xyn[5] = {0, 1 ,2, 3, 4};
    const int xoz[5] = {4, 11, 12, 13, 18};
    const int xnz[5] = {0, 5, 7, 6, 14};
    const int oyz[5] = {3, 7, 10, 13, 17};
    const int nyz[5] = {1, 5, 8, 11, 15};

    //inner cells
    for (i = 1; i <= xlength; i++)
    {
        for (j = 1; j <= xlength; j++)
        {
            for(k = 1; k <= xlength; k++)
            {
                flag = flagField[FINDEXOF(xlength, i, j, k)];
                if (flag == NO_SLIP)
                {
                    for (l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        collideField[INDEXOF(xlength,i, j, k, l)] = FINV(xlength, i, j, k, l);
                    }
                }
                else if (flag == MOVING_WALL)
                {
                    for (l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        ind = INDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j+LATTICEVELOCITIES[l][1], k+LATTICEVELOCITIES[l][2], 0);
                        computeDensity(collideField+ind, &den);
                        collideField[INDEXOF(xlength,i, j, k, l)] = FINV(xlength, i, j, k, l) + WALL(den, l);
                    }
                }
            }
        }
    }

    //xyo and XYN planes
    for(i = 1; i <= xlength; i++)
    {
        for (j = 1; j <= xlength; j++)
        {
            flag = flagField[FINDEXOF(xlength, i, j, 0)];
            if (flag == NO_SLIP) //update 14, 15, 16, 17, 18
            {
                for (l = 0; l<N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, xyo[l])] = FINV(xlength, i, j, k, xyo[l]);
                }
            }
            else if (flag == MOVING_WALL)
            {

                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xyo[l])] = FINV(xlength, i, j, k, xyo[l]) + WALL(den, xyo[l]);
                }
            }

            flag = flagField[FINDEXOF(xlength, i, j, xlength+1)];
            if (flag == NO_SLIP) //update 0, 1, 2, 3, 4
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, xyn[l])] = FINV(xlength, i, j, k, xyn[l]);
                }
            }
            else if (flag == MOVING_WALL)
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xyn[l]][0], j+LATTICEVELOCITIES[xyn[l]][1], k+LATTICEVELOCITIES[xyn[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xyn[l])] = FINV(xlength, i, j, k, xyn[l]) + WALL(den, xyn[l]);
                }
            }
        }
    }

    //OYZ and NYZ
    for(j = 1; j <= xlength; j++)
    {
        for (k = 1; k <= xlength; k++)
        {
            //0YZ
            flag = flagField[FINDEXOF(xlength, 0, j, k)];
            if (flag == NO_SLIP) //update 3, 7, 10, 13, 17
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, oyz[l])] = FINV(xlength, i, j, k, oyz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 3, 7, 10, 13, 17
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, oyz[l])] = FINV(xlength, i, j, k, oyz[l]) + WALL(den, oyz[l]);
                }
            }

            //NYZ
            flag = flagField[FINDEXOF(xlength, xlength+1, j, k)];
            if (flag == NO_SLIP) //update 1, 5, 8, 11, 15
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, nyz[l])] = FINV(xlength, i, j, k, nyz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 1, 5, 8, 11, 15
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[nyz[l]][0], j+LATTICEVELOCITIES[nyz[l]][1], k+LATTICEVELOCITIES[nyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, nyz[l])] = FINV(xlength, i, j, k, nyz[l]) + WALL(den, nyz[l]);
                }
            }
        }
    }

    //X0Z and XNZ
    for(i = 1; i <= xlength; i++)
    {
        for (k = 1; k <= xlength; k++)
        {
            //X0Z
            flag = flagField[FINDEXOF(xlength, i, 0, k)];
            if (flag == NO_SLIP) //update 4, 11, 12, 13, 18
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, xoz[l])] = FINV(xlength, i, j, k, xoz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 4, 11, 12, 13, 18
            {

                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xoz[l]][0], j+LATTICEVELOCITIES[xoz[l]][1], k+LATTICEVELOCITIES[xoz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xoz[l])] = FINV(xlength, i, j, k, xoz[l]) + WALL(den, xoz[l]);
                }
            }

            //XNZ
            flag = flagField[FINDEXOF(xlength, i, xlength+1, k)];
            if (flag == NO_SLIP) //update 0, 5, 7, 6, 14
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, xnz[l])] = FINV(xlength, i, j, k, xnz[l]);
                }
            }
            else if (flag == MOVING_WALL) //update directions 0, 5, 7, 6, 14
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xnz[l]][0], j+LATTICEVELOCITIES[xnz[l]][1], k+LATTICEVELOCITIES[xnz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xnz[l])] = FINV(xlength, i, j, k, xnz[l]) + WALL(den, xnz[l]);
                }
            }
        }
    }


        /* ----- */
        /* Edges */
        /* ----- */

    //X-edges
    for(i = 1; i <= xlength; i++)
    {
        //(i, 0, 0)
        flag = flagField[FINDEXOF(xlength, i, 0, 0)];
        if (flag == NO_SLIP) //update 18
        {
            collideField[INDEXOF(xlength,i, 0, 0, 18)] = FINV(xlength, i, 0, 0, 18);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[18][0], LATTICEVELOCITIES[18][1], LATTICEVELOCITIES[18][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, 0, 0, 18)] = FINV(xlength, i, 0, 0, 18) + WALL(den, 18);
        }

        //(i, 0, xlength+1)
        flag = flagField[FINDEXOF(xlength, i, 0, xlength+1)];
        if (flag == NO_SLIP) //update 4
        {
            collideField[INDEXOF(xlength,i, 0, xlength+1, 4)] = FINV(xlength, i, 0, xlength+1, 4);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[4][0], LATTICEVELOCITIES[4][1], LATTICEVELOCITIES[4][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, 0, xlength+1, 4)] = FINV(xlength, i, 0, xlength+1, 4) + WALL(den, 4);
        }

        //(i, xlength+1, 0)
        flag = flagField[FINDEXOF(xlength, i, xlength+1, 0)];
        if (flag == NO_SLIP) //update 14
        {
            collideField[INDEXOF(xlength,i, xlength+1, 0, 14)] = FINV(xlength, i, xlength+1, 0, 14);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[14][0], xlength+1+LATTICEVELOCITIES[14][1], LATTICEVELOCITIES[14][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, xlength+1, 0, 14)] = FINV(xlength, i, xlength+1, 0, 14) + WALL(den, 14);
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
            collideField[INDEXOF(xlength, i, xlength+1, xlength+1, 0)] = FINV(xlength, i, xlength+1, xlength+1, 0) + WALL(den, 0);
        }
    }


    //Y-edges
    for(j = 1; j <= xlength; j++)
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
            collideField[INDEXOF(xlength, 0, j, 0, 17)] = FINV(xlength, 0, j, 0, 17) + WALL(den, 17);
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
            collideField[INDEXOF(xlength, xlength+1, j, 0, 15)] = FINV(xlength, xlength+1, j, 0, 15) + WALL(den, 15);
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
            collideField[INDEXOF(xlength, 0, j, xlength+1, 3)] = FINV(xlength, 0, j, xlength+1, 3) + WALL(den, 3);
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
            collideField[INDEXOF(xlength,xlength+1, j, xlength+1, 1)] = FINV(xlength, xlength+1, j, xlength+1, 1) + WALL(den, 1);
        }
    }


    //Z-edges
    for(k = 1; k <= xlength; k++)
    {
        //(0, 0, k)
        flag = flagField[FINDEXOF(xlength, 0, 0, k)];
        if (flag == NO_SLIP) //update 13
        {
            collideField[INDEXOF(xlength,0, 0, k, 13)] = FINV(xlength, 0, 0, k, 13);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[13][0], LATTICEVELOCITIES[13][1], k+LATTICEVELOCITIES[13][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, 0, k, 13)] = FINV(xlength, 0, 0, k, 13) + WALL(den, 13);
        }

        //(xlength+1, 0, k)
        flag = flagField[FINDEXOF(xlength, xlength+1, j, 0)];
        if (flag == NO_SLIP) //update 11
        {
            collideField[INDEXOF(xlength,xlength+1, 0, k, 11)] = FINV(xlength, xlength+1, 0, k, 11);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[11][0], LATTICEVELOCITIES[11][1], k+LATTICEVELOCITIES[11][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, 0, k, 11)] = FINV(xlength, xlength+1, 0, k, 11) + WALL(den, 11);
        }

        //(0, xlength+1, k)
        flag = flagField[FINDEXOF(xlength, 0, j, xlength+1)];
        if (flag == NO_SLIP) //update 7
        {
            collideField[INDEXOF(xlength,0, xlength+1, k, 7)] = FINV(xlength, 0, xlength+1, k, 7);
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[7][0], xlength+1+LATTICEVELOCITIES[7][1], k+LATTICEVELOCITIES[7][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, xlength+1, k, 7)] = FINV(xlength, 0, xlength+1, k, 7) + WALL(den, 7);
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
            collideField[INDEXOF(xlength,xlength+1, xlength+1, k, 5)] = FINV(xlength, xlength+1, xlength+1, k, 5) + WALL(den, 5);
        }
    }
    
}

