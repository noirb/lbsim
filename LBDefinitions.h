#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

<<<<<<< HEAD
#define NUMBER_OF_LATTICE_DIRECTIONS 19
#define NUMBER_OF_COORDINATES 3

=======
>>>>>>> 6945c0893aa11e4a2a6c4a26b6cabd61d23f01cb
  static const int LATTICEVELOCITIES[19][3] = {
                    { 0, -1, -1},
                    {-1,  0, -1},
                    { 0,  0, -1},
                    { 1,  0, -1},
                    { 0,  1, -1},
                    {-1, -1,  0},
                    { 0, -1,  0},
                    { 1, -1,  0},
                    {-1,  0,  0},
                    { 0,  0,  0},
                    { 1,  0,  0},
                    {-1,  1,  0},
                    { 0,  1,  0},
                    { 1,  1,  0},
                    { 0, -1,  1},
                    {-1,  0,  1},
                    { 0,  0,  1},
                    { 1,  0,  1},
                    { 0,  1,  1}
  };

  static const double LATTICEWEIGHTS[19] = {
                     1.0 / 36.0,
                     1.0 / 36.0,
                     2.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0,
                     2.0 / 36.0,
                     1.0 / 36.0,
                     2.0 / 36.0,
                    12.0 / 36.0,
                     2.0 / 36.0,
                     1.0 / 36.0,
                     2.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0,
                     2.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0
  };
<<<<<<< HEAD
  static const double C_S = 0.577350269; // 1 / sqrt(3.0);
=======
  static const double C_S = 1.0 / 1.73205080757; // 1 / sqrt(3.0);
>>>>>>> 6945c0893aa11e4a2a6c4a26b6cabd61d23f01cb

  typedef enum { FLUID, NO_SLIP, MOVING_WALL } cell_flag;


  /* ---------------------------------------- */
  /* macros simplify indexing into our arrays */
  /* ---------------------------------------- */

  // used for collide & stream fields to get a distribution at x,y,z,i
  #define INDEXOF(xlength, x, y, z, i) 19 * (z * xlength * xlength + y * xlength + x) + i

  // used for flag field to get flag at x,y,z
  #define FINDEXOF(xlength, x, y, z) (z * xlength * xlength + y * xlength + x)

#endif

