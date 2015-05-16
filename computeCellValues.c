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
    double u_dot_u = dot(velocity, velocity);

    for (int i = 0; i < NUMBER_OF_LATTICE_DIRECTIONS; i++){
        double wp = LATTICEWEIGHTS[i] * *density;
        double u_dot_ci = dot2(LATTICEVELOCITIES[i], velocity);
        feq[i] = 1;

        feq[i] += u_dot_ci / (C_S * C_S);
        feq[i] += (u_dot_ci * u_dot_ci) / (2 * C_S * C_S * C_S * C_S);
        feq[i] += u_dot_u / (2 * C_S * C_S);

        feq[i] *= wp;
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
