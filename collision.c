#include "collision.h"


void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    for( int i = 0; i < NUMBER_OF_LATTICE_DIRECTIONS; i++){
        *(currentCell+i) -= (1.0/(*tau)) * (*(currentCell+i) - *(feq+i));
    }
}

void doCollision(double *collideField, int *flagField, const double * const tau, int xlength){
    double density;
    double velocity[3];
    double feq[19];

    for (int i = 1; i <= xlength; i++){
        for (int j = 1; j <= xlength; j++){
            for (int k = 1; k <= xlength; k++){
                computeDensity(&(collideField[INDEXOF(xlength, i, j, k, 0)]), &density);
                computeVelocity(&(collideField[INDEXOF(xlength, i, j, k, 0)]), &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(&(collideField[INDEXOF(xlength, i, j, k, 0)]), tau, feq);
            }
        }
    }
}
