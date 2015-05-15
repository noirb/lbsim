#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){
    double *collideField=NULL;
    double *streamField=NULL;
    int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    int timesteps;
    int timestepsPerPlotting;

    if (readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv) == -1)
    {
        return 0; // exit if readParameters returned an error
    }

    // collideField & streamField each have Q * (xlength + 2) ^ D cells
    //                    Q == 19, D == 3
    collideField = (double*) malloc(sizeof(double) * 19 * pow(xlength+2, 3));
    streamField  = (double*) malloc(sizeof(double) * 19 * pow(xlength+2, 3));
      // flagField contains (xlength + 2) ^ D cells
    flagField = (int*) malloc(sizeof(int) * pow(xlength+2, 3));
    initialiseFields(collideField, streamField, flagField, xlength);

    for(int t = 0; t < timesteps; t++)
    {
        double *swap=NULL;
        doStreaming(collideField, streamField, flagField, xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField, flagField, &tau, xlength);
        treatBoundary(collideField,flagField,velocityWall,xlength);

        if (t%timestepsPerPlotting==0)
        {
            writeVtkOutput(collideField, flagField, argv[1], t, xlength);
        }
    }

  // cleanup
  free(collideField);
  free(streamField);
  free(flagField);

  return 0;
}

#endif

