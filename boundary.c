#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

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
                        collideField[INDEXOF(xlength,i, j, k, l)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j+LATTICEVELOCITIES[l][1], k+LATTICEVELOCITIES[l][2], 18-l)];
                    }
                }
                else if (flag == MOVING_WALL)
                {
                    for (l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
                    {
                        ind = INDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j+LATTICEVELOCITIES[l][1], k+LATTICEVELOCITIES[l][2], 0);
                        computeDensity(collideField+ind, &den);
                        collideField[INDEXOF(xlength,i, j, k, l)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[l][0], j+LATTICEVELOCITIES[l][1], k+LATTICEVELOCITIES[l][2], 18-l)] + 2*LATTICEWEIGHTS[l]*den*(LATTICEVELOCITIES[l][0]*wallVelocity[0]+LATTICEVELOCITIES[l][1]*wallVelocity[1]+LATTICEVELOCITIES[l][2]*wallVelocity[2])/(C_S*C_S);
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
                    collideField[INDEXOF(xlength,i, j, 0, xyo[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xyo[l]][0], j+LATTICEVELOCITIES[xyo[l]][1], k+LATTICEVELOCITIES[xyo[l]][2], 18-xyo[l])];
                }
            }
            else if (flag == MOVING_WALL)
            {

                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xyo[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xyo[l]][0], j+LATTICEVELOCITIES[xyo[l]][1], k+LATTICEVELOCITIES[xyo[l]][2], 18-xyo[l])] + 2*LATTICEWEIGHTS[xyo[l]]*den*(LATTICEVELOCITIES[xyo[l]][0]*wallVelocity[0]+LATTICEVELOCITIES[xyo[l]][1]*wallVelocity[1]+LATTICEVELOCITIES[xyo[l]][2]*wallVelocity[2])/(C_S*C_S);
                }
            }

            flag = flagField[FINDEXOF(xlength, i, j, xlength+1)];
            if (flag == NO_SLIP) //update 0, 1, 2, 3, 4
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, xyn[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xyn[l]][0], j+LATTICEVELOCITIES[xyn[l]][1], k+LATTICEVELOCITIES[xyn[l]][2], 18-xyn[l])];
                }
            }

            if (flag == MOVING_WALL)
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xyn[l]][0], j+LATTICEVELOCITIES[xyn[l]][1], k+LATTICEVELOCITIES[xyn[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xyn[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xyn[l]][0], j+LATTICEVELOCITIES[xyn[l]][1], k+LATTICEVELOCITIES[xyn[l]][2], 18-xyn[l])] + 2*LATTICEWEIGHTS[xyn[l]]*den*(LATTICEVELOCITIES[xyn[l]][0]*wallVelocity[0]+LATTICEVELOCITIES[xyn[l]][1]*wallVelocity[1]+LATTICEVELOCITIES[xyn[l]][2]*wallVelocity[2])/(C_S*C_S);
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
                for (l = 0; l<N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, oyz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 18-oyz[l])];
                }
            }
            else if (flag == MOVING_WALL) //update directions 3, 7, 10, 13, 17
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, oyz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[oyz[l]][0], j+LATTICEVELOCITIES[oyz[l]][1], k+LATTICEVELOCITIES[oyz[l]][2], 18-oyz[l])] + 2*LATTICEWEIGHTS[oyz[l]]*den*(LATTICEVELOCITIES[oyz[l]][0]*wallVelocity[0]+LATTICEVELOCITIES[oyz[l]][1]*wallVelocity[1]+LATTICEVELOCITIES[oyz[l]][2]*wallVelocity[2])/(C_S*C_S);
                }
            }

            //NYZ
            flag = flagField[FINDEXOF(xlength, xlength+1, j, k)];
            if (flag == NO_SLIP) //update 1, 5, 8, 11, 15
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, nyz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[nyz[l]][0], j+LATTICEVELOCITIES[nyz[l]][1], k+LATTICEVELOCITIES[nyz[l]][2], 18-nyz[l])];
                }
            }
            else if (flag == MOVING_WALL) //update directions 1, 5, 8, 11, 15
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[nyz[l]][0], j+LATTICEVELOCITIES[nyz[l]][1], k+LATTICEVELOCITIES[nyz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, nyz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[nyz[l]][0], j+LATTICEVELOCITIES[nyz[l]][1], k+LATTICEVELOCITIES[nyz[l]][2], 18-nyz[l])] + 2*LATTICEWEIGHTS[nyz[l]]*den*(LATTICEVELOCITIES[nyz[l]][0]*wallVelocity[0]+LATTICEVELOCITIES[nyz[l]][1]*wallVelocity[1]+LATTICEVELOCITIES[nyz[l]][2]*wallVelocity[2])/(C_S*C_S);
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
                    collideField[INDEXOF(xlength,i, j, 0, xoz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xoz[l]][0], j+LATTICEVELOCITIES[xoz[l]][1], k+LATTICEVELOCITIES[xoz[l]][2], 18-xoz[l])];
                }
            }
            else if (flag == MOVING_WALL) //update directions 4, 11, 12, 13, 18
            {

                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xoz[l]][0], j+LATTICEVELOCITIES[xoz[l]][1], k+LATTICEVELOCITIES[xoz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xoz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xoz[l]][0], j+LATTICEVELOCITIES[xoz[l]][1], k+LATTICEVELOCITIES[xoz[l]][2], 18-xoz[l])] + 2*LATTICEWEIGHTS[xoz[l]]*den*(LATTICEVELOCITIES[xoz[l]][0]*wallVelocity[0]+LATTICEVELOCITIES[xoz[l]][1]*wallVelocity[1]+LATTICEVELOCITIES[xoz[l]][2]*wallVelocity[2])/(C_S*C_S);
                }
            }

            //XNZ
            flag = flagField[FINDEXOF(xlength, i, xlength+1, k)];
            if (flag == NO_SLIP) //update 0, 5, 7, 6, 14
            {
                for (l = 0; l < N; l++)
                {
                    collideField[INDEXOF(xlength,i, j, 0, xnz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xnz[l]][0], j+LATTICEVELOCITIES[xnz[l]][1], k+LATTICEVELOCITIES[xnz[l]][2], 18-xnz[l])];
                }
            }
            else if (flag == MOVING_WALL) //update directions 0, 5, 7, 6, 14
            {
                for(l = 0; l < N; l++)
                {
                    ind = INDEXOF(xlength, i + LATTICEVELOCITIES[xnz[l]][0], j+LATTICEVELOCITIES[xnz[l]][1], k+LATTICEVELOCITIES[xnz[l]][2], 0);
                    computeDensity(collideField+ind, &den);
                    collideField[INDEXOF(xlength,i, j, 0, xnz[l])] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[xnz[l]][0], j+LATTICEVELOCITIES[xnz[l]][1], k+LATTICEVELOCITIES[xnz[l]][2], 18-xnz[l])] + 2*LATTICEWEIGHTS[xnz[l]]*den*(LATTICEVELOCITIES[xnz[l]][0]*wallVelocity[0]+LATTICEVELOCITIES[xnz[l]][1]*wallVelocity[1]+LATTICEVELOCITIES[xnz[l]][2]*wallVelocity[2])/(C_S*C_S);
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
            collideField[INDEXOF(xlength,i, 0, 0, 18)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[18][0], LATTICEVELOCITIES[18][1], LATTICEVELOCITIES[18][2], 0)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[18][0], LATTICEVELOCITIES[18][1], LATTICEVELOCITIES[18][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, 0, 0, 18)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[18][0], LATTICEVELOCITIES[18][1], LATTICEVELOCITIES[18][2], 0)] + 2*LATTICEWEIGHTS[18]*den*(LATTICEVELOCITIES[18][0]*wallVelocity[0]+LATTICEVELOCITIES[18][1]*wallVelocity[1]+LATTICEVELOCITIES[18][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(i, 0, xlength+1)
        flag = flagField[FINDEXOF(xlength, i, 0, xlength+1)];
        if (flag == NO_SLIP) //update 4
        {
            collideField[INDEXOF(xlength,i, 0, xlength+1, 4)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[4][0], LATTICEVELOCITIES[4][1], xlength+1+ LATTICEVELOCITIES[4][2], 14)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[4][0], LATTICEVELOCITIES[4][1], LATTICEVELOCITIES[4][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, j, 0, 4)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[4][0], LATTICEVELOCITIES[4][1], xlength+1+LATTICEVELOCITIES[4][2], 14)] + 2*LATTICEWEIGHTS[4]*den*(LATTICEVELOCITIES[4][0]*wallVelocity[0]+LATTICEVELOCITIES[4][1]*wallVelocity[1]+LATTICEVELOCITIES[4][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(i, xlength+1, 0)
        flag = flagField[FINDEXOF(xlength, i, xlength+1, 0)];
        if (flag == NO_SLIP) //update 14
        {
            collideField[INDEXOF(xlength,i, xlength+1, 0, 14)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[14][0], xlength+1+LATTICEVELOCITIES[14][1], LATTICEVELOCITIES[14][2], 4)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[14][0], xlength+1+LATTICEVELOCITIES[14][1], LATTICEVELOCITIES[14][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, xlength+1, 0, 14)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[14][0], xlength+1+LATTICEVELOCITIES[14][1], LATTICEVELOCITIES[14][2], 4)] + 2*LATTICEWEIGHTS[14]*den*(LATTICEVELOCITIES[14][0]*wallVelocity[0]+LATTICEVELOCITIES[14][1]*wallVelocity[1]+LATTICEVELOCITIES[14][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(i, xlength+1, xlength+1)
        flag = flagField[FINDEXOF(xlength, i, xlength+1, xlength+1)];
        if (flag == NO_SLIP) //update 0
        {
            collideField[INDEXOF(xlength,i, xlength+1, xlength+1, 0)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[0][0], xlength+1+LATTICEVELOCITIES[0][1], xlength+1+LATTICEVELOCITIES[0][2], 18)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, i + LATTICEVELOCITIES[0][0], xlength+1+LATTICEVELOCITIES[0][1], xlength+1+LATTICEVELOCITIES[0][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,i, xlength+1, xlength+1, 0)] = collideField[INDEXOF(xlength, i + LATTICEVELOCITIES[0][0], xlength+1+LATTICEVELOCITIES[0][1], xlength+1+LATTICEVELOCITIES[0][2], 18)] + 2*LATTICEWEIGHTS[0]*den*(LATTICEVELOCITIES[0][0]*wallVelocity[0]+LATTICEVELOCITIES[0][1]*wallVelocity[1]+LATTICEVELOCITIES[0][2]*wallVelocity[2])/(C_S*C_S);
        }
    }


    //Y-edges
    for(j = 1; j <= xlength; j++)
    {
        //(0, j, 0)
        flag = flagField[FINDEXOF(xlength, 0, j, 0)];
        if (flag == NO_SLIP) //update 17
        {
            collideField[INDEXOF(xlength,0, j, 0, 17)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[17][0], j+LATTICEVELOCITIES[17][1], LATTICEVELOCITIES[17][2], 1)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[17][0], j+LATTICEVELOCITIES[17][1], LATTICEVELOCITIES[17][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, j, 0, 17)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[17][0], j+LATTICEVELOCITIES[17][1], LATTICEVELOCITIES[17][2], 1)] + 2*LATTICEWEIGHTS[17]*den*(LATTICEVELOCITIES[17][0]*wallVelocity[0]+LATTICEVELOCITIES[17][1]*wallVelocity[1]+LATTICEVELOCITIES[17][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(xlength+1, j, 0)
        flag = flagField[FINDEXOF(xlength, xlength+1, j, 0)];
        if (flag == NO_SLIP) //update 15
        {
            collideField[INDEXOF(xlength,xlength+1, j, 0, 15)] = collideField[INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[15][0], j+LATTICEVELOCITIES[15][1], LATTICEVELOCITIES[15][2], 3)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[15][0], j+LATTICEVELOCITIES[15][1], LATTICEVELOCITIES[15][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, j, 0, 15)] = collideField[INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[15][0],j +LATTICEVELOCITIES[15][1], LATTICEVELOCITIES[15][2], 3)] + 2*LATTICEWEIGHTS[15]*den*(LATTICEVELOCITIES[15][0]*wallVelocity[0]+LATTICEVELOCITIES[15][1]*wallVelocity[1]+LATTICEVELOCITIES[15][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(0, j, xlength+1)
        flag = flagField[FINDEXOF(xlength, 0, j, xlength+1)];
        if (flag == NO_SLIP) //update 3
        {
            collideField[INDEXOF(xlength,0, j, xlength+1, 3)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[3][0], j+LATTICEVELOCITIES[3][1], xlength+1+LATTICEVELOCITIES[3][2], 15)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[3][0], j+LATTICEVELOCITIES[3][1], xlength+1+LATTICEVELOCITIES[3][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, j, xlength+1, 3)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[3][0], j+LATTICEVELOCITIES[3][1], xlength+1+LATTICEVELOCITIES[3][2], 15)] + 2*LATTICEWEIGHTS[3]*den*(LATTICEVELOCITIES[3][0]*wallVelocity[0]+LATTICEVELOCITIES[3][1]*wallVelocity[1]+LATTICEVELOCITIES[3][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(xlength+1, j, xlength+1)
        flag = flagField[FINDEXOF(xlength, xlength+1, j, xlength+1)];
        if (flag == NO_SLIP) //update 1
        {
            collideField[INDEXOF(xlength,xlength+1, j, xlength+1, 1)] = collideField[INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[1][0], j+LATTICEVELOCITIES[1][1], xlength+1+LATTICEVELOCITIES[1][2], 17)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[1][0], j+LATTICEVELOCITIES[1][1], xlength+1+LATTICEVELOCITIES[1][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, j, xlength+1, 1)] = collideField[INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[1][0], j+LATTICEVELOCITIES[1][1], xlength+1+LATTICEVELOCITIES[1][2], 17)] + 2*LATTICEWEIGHTS[0]*den*(LATTICEVELOCITIES[1][0]*wallVelocity[0]+LATTICEVELOCITIES[1][1]*wallVelocity[1]+LATTICEVELOCITIES[1][2]*wallVelocity[2])/(C_S*C_S);
        }
    }


    //Z-edges
    for(k = 1; k <= xlength; k++)
    {
        //(0, 0, k)
        flag = flagField[FINDEXOF(xlength, 0, 0, k)];
        if (flag == NO_SLIP) //update 13
        {
            collideField[INDEXOF(xlength,0, 0, k, 13)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[13][0], LATTICEVELOCITIES[13][1], k+LATTICEVELOCITIES[13][2], 5)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[13][0], LATTICEVELOCITIES[13][1], k+LATTICEVELOCITIES[13][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, 0, k, 13)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[13][0], LATTICEVELOCITIES[13][1], k+LATTICEVELOCITIES[13][2], 5)] + 2*LATTICEWEIGHTS[13]*den*(LATTICEVELOCITIES[13][0]*wallVelocity[0]+LATTICEVELOCITIES[13][1]*wallVelocity[1]+LATTICEVELOCITIES[13][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(xlength+1, 0, k)
        flag = flagField[FINDEXOF(xlength, xlength+1, j, 0)];
        if (flag == NO_SLIP) //update 11
        {
            collideField[INDEXOF(xlength,xlength+1, 0, k, 11)] = collideField[INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[11][0], LATTICEVELOCITIES[11][1], k+LATTICEVELOCITIES[11][2], 7)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[11][0], LATTICEVELOCITIES[11][1], k+LATTICEVELOCITIES[11][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, 0, k, 11)] = collideField[INDEXOF(xlength, xlength+1 + LATTICEVELOCITIES[11][0],LATTICEVELOCITIES[11][1], k+LATTICEVELOCITIES[11][2], 7)] + 2*LATTICEWEIGHTS[11]*den*(LATTICEVELOCITIES[11][0]*wallVelocity[0]+LATTICEVELOCITIES[11][1]*wallVelocity[1]+LATTICEVELOCITIES[11][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(0, xlength+1, k)
        flag = flagField[FINDEXOF(xlength, 0, j, xlength+1)];
        if (flag == NO_SLIP) //update 7
        {
            collideField[INDEXOF(xlength,0, xlength+1, k, 7)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[7][0], xlength+1+LATTICEVELOCITIES[7][1], k+LATTICEVELOCITIES[7][2], 11)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, LATTICEVELOCITIES[7][0], xlength+1+LATTICEVELOCITIES[7][1], k+LATTICEVELOCITIES[7][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,0, xlength+1, k, 7)] = collideField[INDEXOF(xlength, LATTICEVELOCITIES[7][0], xlength+1+LATTICEVELOCITIES[7][1], k+LATTICEVELOCITIES[7][2], 11)] + 2*LATTICEWEIGHTS[7]*den*(LATTICEVELOCITIES[7][0]*wallVelocity[0]+LATTICEVELOCITIES[7][1]*wallVelocity[1]+LATTICEVELOCITIES[7][2]*wallVelocity[2])/(C_S*C_S);
        }

        //(xlength+1, xlength+1, k)
        flag = flagField[FINDEXOF(xlength, xlength+1, xlength+1, k)];
        if (flag == NO_SLIP) //update 5
        {
            collideField[INDEXOF(xlength,xlength+1, xlength+1, k, 5)] = collideField[INDEXOF(xlength, xlength+1+LATTICEVELOCITIES[5][0], xlength+1+LATTICEVELOCITIES[5][1], k+LATTICEVELOCITIES[5][2], 11)];
        }
        else if (flag == MOVING_WALL)
        {
            ind = INDEXOF(xlength, xlength+1+LATTICEVELOCITIES[5][0], xlength+1+LATTICEVELOCITIES[5][1], k+LATTICEVELOCITIES[5][2], 0);
            computeDensity(collideField+ind, &den);
            collideField[INDEXOF(xlength,xlength+1, xlength+1, k, 5)] = collideField[INDEXOF(xlength, xlength+1+LATTICEVELOCITIES[5][0], xlength+1+LATTICEVELOCITIES[5][1], k+LATTICEVELOCITIES[5][2], 11)] + 2*LATTICEWEIGHTS[5]*den*(LATTICEVELOCITIES[5][0]*wallVelocity[0]+LATTICEVELOCITIES[5][1]*wallVelocity[1]+LATTICEVELOCITIES[5][2]*wallVelocity[2])/(C_S*C_S);
        }
    }
    
}

