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
    
    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);
    // TODO: initialise pointers here!
    
    initialiseFields(collideField,streamField,flagField,xlength);


  return 0;
}

#endif

