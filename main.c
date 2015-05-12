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

    initialiseFields(collideField, streamField, flagField, xlength);

/*    
    for(int t = 0; t < timesteps; t++)
    {
        double *swap=NULL;
        doStreaming(collideField,streamField,flagfield,xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        doCollision(collideField,flagfield,&tau,xlength);
        treatBoundary(collideField,flagfield,velocityWall,xlength);
        
        if (t%timestepsPerPlotting==0)
        {
            writeVtkOutput(collideField,flagfield,argv,t,xlength);
        }
    }
*/

  return 0;
}

#endif

