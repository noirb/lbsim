#ifndef _MAIN_C_
#define _MAIN_C_

#ifdef _WIN32
#define NOMINMAX
#endif

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){
    double    *collideField=NULL;
    double    *streamField=NULL;
    flag_data *flagField=NULL;
    int xlength;
    int ylength;
    int zlength;
    int totalElements;
    double tau;
    int timesteps;
    int timestepsPerPlotting;

    if (readParameters(&xlength, &ylength, &zlength, &tau, &timesteps, &timestepsPerPlotting, argc, argv) == -1)
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
    flagField = (flag_data*) calloc(totalElements, sizeof(flag_data));

    // exit if any allocation failed
    if (collideField == NULL || streamField == NULL || flagField == NULL)
    {
        fprintf(stderr, "Could not allocate space for %d elements!\n", totalElements);
        return 0;
    }

    printf("Initializing fields..."); fflush(stdout);
    if (initialiseFields(collideField, streamField, flagField, xlength, ylength, zlength, argv[1]) == -1)
    {
        fprintf(stderr, "Errors encountered in input. Could not initialize domain. Exiting...\n");
        return 0;
    }
    printf("Done!\n");

    for(int t = 0; t < timesteps; t++)
    {
        double *swap=NULL;
        doStreaming(collideField, streamField, flagField, xlength, ylength, zlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;

        doCollision(collideField, flagField, &tau, xlength, ylength, zlength);
        treatBoundary(collideField, flagField, xlength, ylength, zlength);

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

