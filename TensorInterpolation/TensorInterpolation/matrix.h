#ifndef _MATRIX_
#define _MATRIX_
#define _USE_MATH_DEFINES

#include <iostream>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include "exp-map.h"

/*matrix mul a*b=c */
void gsl_matrix_mul(gsl_matrix *a, gsl_matrix *b, gsl_matrix *c);
/* norm */
double normMatrix( const gsl_matrix* src );
/* trace */
double traceMatrix( const gsl_matrix* src );
/* normalize */
void normalizeMatrix( gsl_matrix* dest, const gsl_matrix* src );
/* inner product */
double innerProduct( const gsl_matrix* m1, const gsl_matrix* m2 );
/* inverse matrix Äæ¾ØÕó*/
void invertMatrix( gsl_matrix* dest, const gsl_matrix* src );
/* square root of matrix */
void sqrtMatrix( gsl_matrix* dest, const gsl_matrix* src );
/* logarithm of matrix */
void logMatrix( gsl_matrix* dest, const gsl_matrix* src );
/* exponential of matrix */
void expMatrix( gsl_matrix* dest, const gsl_matrix* src );
/* interpolate matrix */
void interpolateMatrix( gsl_matrix* dest, const gsl_matrix* src, double param );
/* logarithm of symmetrical matrix */
void logMatrix_sym( gsl_matrix* dest, const gsl_matrix* src );
/* exponential of symmetrical matrix */
void expMatrix_sym( gsl_matrix* dest, const gsl_matrix* src );
/* isotropic of symmetrical matrix */
void isotropicMatrix_sym( gsl_matrix* dest, const gsl_matrix* src );
/* deviatoric of symmetrical matrix */
void deviatoricMatrix_sym( gsl_matrix* dest, const gsl_matrix* src );
void eigenSolve_sym( gsl_vector* eigVal, gsl_matrix* eigVec, const gsl_matrix* src, const bool zeroCheck );
void invEigenSolve_sym( gsl_matrix* dest, const gsl_vector* eigVal, const gsl_matrix* eigVec );
double detMatrix( const gsl_matrix* src );
void rotationMatrix( gsl_matrix* dest, const gsl_vector* axis, const double degree );
void getRotMatrix( gsl_matrix* dest, const gsl_matrix* src1, const gsl_matrix* src2 );
void interpolateRotMatrix( gsl_matrix* dest, const gsl_matrix* src, const double ratio );
void expMap( gsl_matrix* rotMatrix, const gsl_vector* axis, const double angle );

#endif
