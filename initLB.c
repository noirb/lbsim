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
  for (int i = 0; i < xlength+2; i++)
  {
    for (int j = 0; j < xlength+2; j++)
    {
      for (int k = 0; k < xlength+2; k++)
      {
        // set distributions at (i,j,k)
        if (i == 0 || i == xlength+1 || j == 0 || j == xlength+1 || k == 0 || k == xlength+1){
            for (int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
            {
                collideField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
                streamField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
                //collideField[INDEXOF(xlength, i, j, k, l)] = 0;
                //streamField[INDEXOF(xlength, i, j, k, l)] = 0;
            }
            //collideField[INDEXOF(xlength, i, j, k, 9)] = 1;
            //streamField[INDEXOF(xlength, i, j, k, 9)] = 1;
            flagField[FINDEXOF(xlength, i, j, k)] = NO_SLIP;
        } else {
            for (int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
            {
                collideField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
                streamField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
            }
            flagField[FINDEXOF(xlength, i, j, k)] = FLUID;
        }
      }
    }
  }

  for (int i = 0; i < xlength+2; i++){
    for (int k = 0; k < xlength+2; k++){
        flagField[FINDEXOF(xlength, i, xlength+1, k)] = MOVING_WALL;
    }
  }
}
