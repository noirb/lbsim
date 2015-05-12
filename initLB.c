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
  /* TODO */
}

