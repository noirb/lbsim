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
    char cellDataFile[20];
    int xlength;
    int ylength;
    int zlength;
    int totalElements;
    double tau;
    double velocityWall[3];
    int timesteps;
    int timestepsPerPlotting;

    if (readParameters(&xlength, &ylength, &zlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, cellDataFile, argc, argv) == -1)
    {
        fprintf(stderr, "Error reading parameter file: %s\n", argv[1]);
        return 0; // exit if readParameters returned an error
    }
    totalElements = (xlength+2) * (ylength+2) * (zlength+2);

    // collideField & streamField each have Q * (totalElements) ^ D cells
    //                    Q == 19, D == 3
    collideField = (double*) calloc(19 * totalElements, sizeof(double));
    streamField  = (double*) calloc(19 * totalElements, sizeof(double));
      // flagField contains (totalElements) ^ D cells
    flagField = (int*) calloc(totalElements, sizeof(int));

    // exit if any allocation failed
    if (collideField == NULL || streamField == NULL || flagField == NULL)
    {
        fprintf(stderr, "Could not allocate space for %d elements!\n", totalElements);
        return 0;
    }

    printf("Initializing fields..."); fflush(stdout);
    initialiseFields(collideField, streamField, flagField, xlength, ylength, zlength, cellDataFile);
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
            
            writeVtkOutput(collideField, flagField, argv[1], t, xlength, ylength, zlength);
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

