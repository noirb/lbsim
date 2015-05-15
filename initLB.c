#include "initLB.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
  if (argc != 2)
  {
    fprintf(stderr, "ERROR: Incorrect number of inputs given!\n\tUsage:\n\t %s inputFile", argv[0]);
    return -1;
  }

  READ_INT(argv[1],     *xlength);
  READ_DOUBLE(argv[1],  *tau);
  READ_DOUBLE(argv[1],  *velocityWall);
  READ_INT(argv[1],     *timesteps);
  READ_INT(argv[1],     *timestepsPerPlotting);
  

  return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

        /* Initialize Arrays */
  for (int i = 0; i < xlength; i++)
  {
    for (int j = 0; j < xlength; j++)
    {
      for (int k = 0; k < xlength; k++)
      {
        // set distributions at (i,j,k)
        for (int l = 0; l < xlength; l++)
        {
            collideField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
            streamField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
        }

        // set flags at (i,j,k)
        flagField[FINDEXOF(xlength, i, j, k)] = k == 0 ? NO_SLIP :                      // bottom layer
                                                k == xlength - 1 ? MOVING_WALL :        // top layer
                                                i == 0 || i == xlength - 1 ? NO_SLIP :  // east/west walls
                                                j == 0 || j == xlength - 1 ? NO_SLIP :  // north/south walls
                                                FLUID;                                  // interior cells
      }
    }
  }
}

