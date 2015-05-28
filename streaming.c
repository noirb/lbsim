#include "streaming.h"
#include <stdio.h>

void doStreaming(double *collideField, double *streamField, flag_data *flagField, int xlength, int ylength, int zlength)
{
    for (int i = 1; i <= xlength; i++){
        for (int j = 1; j <= ylength; j++){
            for (int k = 1; k <= zlength; k++){
                if (flagField[FINDEXOF(i, j, k)].flag == FLUID)
                {
                    for (int dir = 0; dir < NUMBER_OF_LATTICE_DIRECTIONS; dir++){
                        streamField[INDEXOF(i, j, k, dir)] = collideField[INDEXOF(i+LATTICEVELOCITIES[18-dir][0], j+LATTICEVELOCITIES[18-dir][1], k+LATTICEVELOCITIES[18-dir][2], dir)];
                    }
                }
            }
        }
    }
}

