#pragma once
#ifndef _MATRIX_
#define _MATRIX_

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

/* logarithm of matrix */
void logMatrix(gsl_matrix* dest, const gsl_matrix* src);
/* exponential of matrix */
void expMatrix(gsl_matrix* dest, const gsl_matrix* src);
/* interpolate matrix */
void interpolateMatrix(gsl_matrix* dest, const gsl_matrix* src, double param);

#endif