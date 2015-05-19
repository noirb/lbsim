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
    collideField = (double*) calloc(19 * pow(xlength+2, 3), sizeof(double));
    streamField  = (double*) calloc(19 * pow(xlength+2, 3), sizeof(double));
      // flagField contains (xlength + 2) ^ D cells
    flagField = (int*) calloc(pow(xlength+2, 3), sizeof(int));

    printf("Initializing fields..."); fflush(stdout);
    initialiseFields(collideField, streamField, flagField, xlength);
    printf("Done!\n");

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
            printf("Writing output for timestep %d/%d...", t, timesteps); fflush(stdout);
            
            writeVtkOutput(collideField, flagField, argv[1], t, xlength);
            printf("Done!\n");
        }
    }

  // cleanup
  free(collideField);
  free(streamField);
  free(flagField);

  return 0;
}

#endif

