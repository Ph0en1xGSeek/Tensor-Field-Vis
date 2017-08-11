#include "matrixLite.h"



void logMatrix(gsl_matrix* dest, const gsl_matrix* src) {
	double threshold;
	int k, i;
	int iterate;

	gsl_matrix* I = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* prevA = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* nextA = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* TMP = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* prevZ = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* nextZ = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* X = gsl_matrix_alloc(src->size1, src->size2);

	// initialize 
	gsl_matrix_set_identity(I);
	gsl_matrix_memcpy(prevA, src);
	gsl_matrix_memcpy(TMP, src);
	gsl_matrix_sub(TMP, I);
	k = 0;
	threshold = normMatrix(TMP);

	iterate = 0;
	while (threshold > 0.5 && iterate < _MAXITERATE) {
		sqrtMatrix(nextA, prevA);
		k++;
		gsl_matrix_swap(prevA, nextA);
		gsl_matrix_memcpy(TMP, prevA);
		gsl_matrix_sub(TMP, I);
		//    if (threshold < normMatrix(TMP))
		//      threshold = 0.4;
		//    else
		threshold = normMatrix(TMP);
		++iterate;
	}
	gsl_matrix_set_identity(nextA);
	gsl_matrix_sub(nextA, prevA);
	gsl_matrix_swap(nextA, prevA);

	gsl_matrix_memcpy(prevZ, prevA);
	gsl_matrix_memcpy(X, prevA);
	i = 1;
	threshold = normMatrix(prevZ);
	iterate = 0;
	while (threshold > _EPS && iterate < _MAXITERATE) {
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, prevZ, prevA, 0.0, nextZ);
		gsl_matrix_swap(nextZ, prevZ);
		gsl_matrix_memcpy(TMP, prevZ);
		i++;
		gsl_matrix_scale(TMP, 1.0 / i);
		gsl_matrix_add(X, TMP);
		//    if(threshold < normMatrix(prevZ))
		//      threshold = 0.0;
		//    else 
		threshold = normMatrix(prevZ);
		++iterate;
	}
	gsl_matrix_scale(X, pow(2, k));

	gsl_matrix_memcpy(dest, X);

	// free memory
	gsl_matrix_free(X);
	gsl_matrix_free(prevZ);
	gsl_matrix_free(nextZ);
	gsl_matrix_free(prevA);
	gsl_matrix_free(nextA);
	gsl_matrix_free(TMP);
	gsl_matrix_free(I);
	return;
}

void expMatrix(gsl_matrix* dest, const gsl_matrix* src) {
	double c;
	int i, j, k;
	const int q = 6; // this value is come from the thesis
	gsl_matrix* A = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* D = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* invD = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* N = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* prevX = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix* nextX = gsl_matrix_alloc(src->size1, src->size2);
	//  gsl_matrix* XX = gsl_matrix_alloc( src->size1, src->size2 );
	gsl_matrix* TMP = gsl_matrix_alloc(src->size1, src->size2);

	// initialize
	gsl_matrix_memcpy(A, src);
	j = 1 + floor(log(normMatrix(A)));
	if (j < 0.0) j = 0.0;
	gsl_matrix_scale(A, pow(2.0, -j));
	gsl_matrix_set_identity(D);
	gsl_matrix_set_identity(N);
	gsl_matrix_set_identity(prevX);
	c = 1.0;

	for (k = 1; k <= q; ++k) {
		c = c * (q - (double)k + 1.0) / ((double)k * (2.0 * q - (double)k + 1.0));
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, prevX, 0.0, nextX);
		gsl_matrix_memcpy(prevX, nextX);
		gsl_matrix_scale(nextX, c);
		gsl_matrix_add(N, nextX);
		gsl_matrix_scale(nextX, pow(-1.0, k));
		gsl_matrix_add(D, nextX);
	}

	invertMatrix(invD, D);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invD, N, 0.0, prevX);
	//  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, prevX, prevX, 0.0, XX );
	gsl_matrix_set_identity(dest);

	for (i = 0; i < (int)pow(2, j); ++i) {
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dest, prevX, 0.0, TMP);
		gsl_matrix_memcpy(dest, TMP);
	}

	gsl_matrix_free(A);
	gsl_matrix_free(D);
	gsl_matrix_free(invD);
	gsl_matrix_free(N);
	gsl_matrix_free(prevX);
	gsl_matrix_free(nextX);
	// gsl_matrix_free( XX );
	gsl_matrix_free(TMP);
	return;
}

void interpolateMatrix(gsl_matrix* dest, const gsl_matrix* src, double ratio) {
	assert(src->size1 == src->size2);
	assert(dest->size1 == src->size1);
	assert(dest->size2 == src->size2);
	if (ratio < 0.0) ratio = 0.0;
	if (ratio > 1.0) ratio = 1.0;

	int i;
	unsigned int bitFlag;

	// make matrix array
	// array[0] = M^(1/2)
	// array[1] = M^(1/4)
	// ...
	// array[n-1] = M^(1/pow(2, n ))
	// const unsigned int intsize = sizeof( int );
	gsl_matrix* tmp = gsl_matrix_alloc(src->size1, src->size2);
	gsl_matrix** roots = new gsl_matrix*[_INTSIZE];
	for (i = 0; i < _INTSIZE; ++i) {
		roots[i] = gsl_matrix_alloc(src->size1, src->size2);
	}
	sqrtMatrix(roots[0], src);
	for (i = 0; i + 1 < _INTSIZE; ++i) {
		sqrtMatrix(roots[i + 1], roots[i]);
	}

	bitFlag = static_cast<unsigned int>(floor(ratio * static_cast<double>(0xffffffff)));
	gsl_matrix_set_identity(dest);
	for (i = _INTSIZE - 1; i >= 0; i--) {
		if ((bitFlag & 0x00000001) != 0) {
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dest, roots[i], 0.0, tmp);
			gsl_matrix_memcpy(dest, tmp);
		}
		bitFlag = bitFlag >> 1;
	}
	for (i = 0; i < _INTSIZE; ++i) {
		gsl_matrix_free(roots[i]);
	}
	delete roots;

	return;
}
