/*
 * gravmag_primitive_subs.h -- header file for primitive vertical prism
 *                             gravmag subroutines
 * 
 *
 * CITATION
 * =======
 *
 * To cite this software in publications, please use:
 *
 *  Cronkite-Ratcliff, C., Phelps, G., Scheirer, D., 2022, gravmagsubs:
 *  Gravitational and magnetic attraction of 3-D vertical rectangular
 *  prisms, U.S. Geological Survey software release, version 1.0,
 *  https://doi.org/10.5066/P9HBSPRM
 *
 *
 */

#include <cmath>
#include <cfloat>

// Sets machine epsilon tolerance
#define eps DBL_EPSILON 

// Sets length of fixed-length char array subname
#define NDIMSUBNAME 50

// universal graviational constant
#define BigGx1e11       6.674          /* units of m**3/kg/sec**2 * 1.E11 */

// radians per degree
#define RAD_PER_DEG   ((M_PI)/180.)

int    rect_prism_grav1(double, double, double, double, double, double, double, double, double, double, double *, int) ;
int    rect_prism_mag1(double, double, double, double, double, double, double, double, double, double, double, double,
                        double, double, double, double *, int);
void   add_tangents(double, double *, double *) ;

