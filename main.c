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

  return 0;
}

#endif

