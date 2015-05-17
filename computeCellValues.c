#include "computeCellValues.h"

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
    double u_dot_u = dot(velocity, velocity);
    double temp = 0;

    for (int i = 0; i < NUMBER_OF_LATTICE_DIRECTIONS; i++){
        feq[i] = 0;
        feq[i] += C_S * C_S * C_S * C_S;
        temp = 0;
        for (int j = 0; j < NUMBER_OF_COORDINATES; j++){
            feq[i] += C_S * C_S * LATTICEVELOCITIES[i][j] * velocity[j];
            temp += LATTICEVELOCITIES[i][j] * velocity[j];
        }
        feq[i] -= C_S * C_S * u_dot_u * 0.5;
        feq[i] += temp * temp * 0.5;
        feq[i] *= LATTICEWEIGHTS[i] * (*density)/(C_S * C_S * C_S * C_S);
    }
}

double dot(const double * const v1, const double * const v2)
{
    double result = 0;
    for (int i = 0; i < NUMBER_OF_COORDINATES; i++)
    {
        result += v1[i] * v2[i];
    }

    return result;
}

double dot2(const int * const v1, const double * const v2)
{
    double result = 0;
    for (int i = 0; i < NUMBER_OF_COORDINATES; i++)
    {
        result += (double)(v1[i]) * v2[i];
    }

    return result;
}
