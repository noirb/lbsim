#include "streaming.h"
#include <stdio.h>

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
int flag;
    for (int i = 1; i <= xlength; i++){
        for (int j = 1; j <= xlength; j++){
            for (int k = 1; k <= xlength; k++){
                for (int dir = 0; dir < NUMBER_OF_LATTICE_DIRECTIONS; dir++){
                    if ((flag = flagField[FINDEXOF(xlength, i, j, k)]) == FLUID)
                        streamField[INDEXOF(xlength, i, j, k, dir)] = collideField[INDEXOF(xlength, i+LATTICEVELOCITIES[18-dir][0], j+LATTICEVELOCITIES[18-dir][1], k+LATTICEVELOCITIES[18-dir][2], dir)];
                }
            }
        }
    }
}

