#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
    for (int i = 1; i <= xlength; i++){
        for (int j = 1; j <= xlength; j++){
            for (int k = 1; k <= xlength; k++){
                for (int dir = 1; dir <= xlength; dir++){
                    streamField[INDEXOF(xlength, i, j, k, dir)] = collideField[INDEXOF(xlength, i+c[dir][0], j+c[dir][1], k+c[dir][2], dir)];
                }
            }
        }
    }
}

