#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

//#define NDEBUG      // uncomment this line to disable asserts in INDEXOF and FINDEXOF. MUST be defined before including assert.h!
#include <assert.h>

// math.h is only used for pow() to compute the max index to check in the asserts below
// don't include it if we're not going to use the asserts, anyway
#ifndef NDEBUG
#include <math.h>
#endif

#define NUMBER_OF_LATTICE_DIRECTIONS 19
#define NUMBER_OF_COORDINATES 3

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

  static const double C_S = 1.0 / 1.73205080757; // 1 / sqrt(3.0);

  typedef enum { FLUID, NO_SLIP, MOVING_WALL } cell_flag;


  /* ---------------------------------------- */
  /* macros simplify indexing into our arrays */
  /* ---------------------------------------- */

  // used for collide & stream fields to get a distribution at x,y,z,i
  #define INDEXOF(xlength, x, y, z, i) __extension__({ \
                                                       int index = _INDEXOF_((xlength), (x), (y), (z), (i)); \
                                                       assert(index >= 0); assert(index < NUMBER_OF_LATTICE_DIRECTIONS * pow(xlength+2, NUMBER_OF_COORDINATES)); \
                                                       index; \
                                                     })
  #define _INDEXOF_(xlength, x, y, z, i) 19 * ((z) * xlength * xlength + (y) * xlength + (x)) + (i)

  // used for flag field to get flag at x,y,z
  #define FINDEXOF(xlength, x, y, z) __extension__({ \
                                                     int index = _FINDEXOF_((xlength), (x), (y), (z)); \
                                                     assert(index >= 0); assert(index < pow(xlength+2, NUMBER_OF_COORDINATES)); \
                                                     index; \
                                                     })
  #define _FINDEXOF_(xlength, x, y, z) ((z) * xlength * xlength + (y) * xlength + (x))

#endif

