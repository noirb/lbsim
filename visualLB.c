#include "helper.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void writeVtkOutput(const double * const collideField, const int * const flagField, const char * filename, unsigned int t, int xlength) {
  int i,j,k;
  char szFileName[80];
  FILE *fp=NULL;
  sprintf( szFileName, "%s.%i.vtk", filename, t );
  fp = fopen( szFileName, "w");

  if( fp == NULL )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
  }

  // allocate space to store densities
  double *densities = (double*) malloc(sizeof(double) * pow(xlength+2, 3));

  write_vtkHeader( fp, xlength);
  write_vtkPointCoordinates(fp, xlength, 1, 1, 1);

  fprintf(fp,"\n");
  fprintf(fp,"CELL_DATA %i \n", (xlength+1)*(xlength+1)*(xlength+1) ); // note: dimensions for cells must be different from points!
  fprintf(fp, "SCALARS density float 1 \n"); 
  fprintf(fp, "LOOKUP_TABLE default \n");
  for(j = 0; j < xlength+2; j++) {
    for(i = 0; i < xlength+2; i++) {
      for(k = 0; k < xlength+2; k++) {
        computeDensity(&collideField[INDEXOF(xlength, i,j,k,0)], &densities[FINDEXOF(xlength, i,j,k)]);
        if (i > 0 && j > 0 && k > 0)
            fprintf(fp, "%f\n", densities[FINDEXOF(xlength, i,j,k)] );
      }
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"POINT_DATA %i \n", (xlength+2)*(xlength+2)*(xlength+2) );
  fprintf(fp, "VECTORS velocity float\n");
  for(j = 0; j < xlength+2; j++) {
    for(i = 0; i < xlength+2; i++) {
      for(k = 0; k < xlength+2; k++) {
        double velocity[3];
        computeVelocity(&collideField[INDEXOF(xlength, i,j,k,0)], densities, velocity);
        fprintf(fp, "%f %f %f\n", velocity[0], velocity[1], velocity[2]);
      }
    }
  }

  free(densities);

  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    ERROR( szBuff );
  }
}

void write_vtkHeader( FILE *fp, int xlength)
{

  if( fp == NULL )
  {
    char szBuff[80];
    sprintf( szBuff, "Null pointer in write_vtkHeader" );
    ERROR( szBuff );
    return;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"\n");	
  fprintf(fp,"DATASET STRUCTURED_GRID\n");
  fprintf(fp,"DIMENSIONS  %i %i %i\n", xlength+2, xlength+2, xlength+2);
  fprintf(fp,"POINTS %i float\n", (xlength+2)*(xlength+2)*(xlength+2) );
  fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int xlength, double dx, double dy, double dz)
{
  double originX = 0.0;
  double originY = 0.0;
  double originZ = 0.0;

  int i = 0;
  int j = 0;
  int k = 0;

  for(i = 0; i < xlength+2; i++) {
    for(j = 0; j < xlength+2; j++) {
      for (k = 0; k < xlength+2; k++) {
        fprintf(fp, "%f %f %f\n", originX+(i*dx), originY+(j*dy), originZ+(k*dz) );
      }
    }
  }
}

