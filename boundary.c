#include "boundary.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  //loop over the boundaries
	int i,j;
	for (i = 0; i<xlength; i++)
		for(j = 0; j<xlength; j++)
	{
		//y-z plane with x = 0
		collidefield(FINDEXOF(0, i, j));

		//
	}
}

