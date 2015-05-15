#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
    *density = 0;
    for (int i = 0; i < NUMBER_OF_LATTICE_DIRECTIONS; i++){
        *density += *(currentCell + i);
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    for (int j = 0; j < NUMBER_OF_COORDINATES; j++){
        velocity[j] = 0;
        for( int i = 0; i < NUMBER_OF_LATTICE_DIRECTIONS; i++){
            velocity[j] += *(currentCell+i) * LATTICEVELOCITIES[i][j];
        }
        velocity[j] /= *density;
    }
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
    double temp, u_dot_u = 0;

    for (int j = 0; j < NUMBER_OF_COORDINATES; j++){
        u_dot_u += velocity[j] * velocity[j];   // Precalculate u.u to save time
    }

    for (int i = 0; i < NUMBER_OF_LATTICE_DIRECTIONS; i++){
        feq[i] = 0;
        feq[i] += C_S * C_S * C_S * C_S;
        temp = 0;
        for (int j = 0; j < NUMBER_OF_COORDINATES; j++){
            feq[i] += C_S * C_S * LATTICEVELOCITIES[i][j] * velocity[j];
            feq[i] -= C_S * C_S * u_dot_u * 0.5;
            temp += LATTICEVELOCITIES[i][j] * velocity[j];
        }
        feq[i] += temp * temp * 0.5;
        feq[i] *= LATTICEWEIGHTS[i] * (*density)/(C_S * C_S * C_S * C_S);
    }
}

