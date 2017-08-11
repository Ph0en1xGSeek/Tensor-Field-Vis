#include "stdafx.h"
#include "matrix.h"
#define M_PI acos(-1)

using namespace std;

#define _EPS      (1.0e-7)
#define _INTSIZE      ( 32)
#define _MAXITERATE  (1000)

void gsl_matrix_mul(gsl_matrix *a, gsl_matrix *b, gsl_matrix *c)
{
	for (size_t i = 0; i<a->size1; i++)
	{
		for (size_t j = 0; j<b->size2; j++)
		{
			double sum = 0.0;
			for (size_t k = 0; k<b->size1; k++)
			{
				sum += gsl_matrix_get(a, i, k)*gsl_matrix_get(b, k, j);
			}
			gsl_matrix_set(c, i, j, sum);
		}
	}
}

double normMatrix( const gsl_matrix* src ){
  double norm = 0.0;
  int i, j;
  for( i = 0; i < (int)src->size1; ++i )
    for( j = 0; j < (int)src->size2; ++j )
      norm += gsl_matrix_get( src, i, j ) * gsl_matrix_get( src, i, j );
  if( norm == 0.0 ){
    return 0.0;
  }
  norm = sqrt( norm );
  return norm;
}

double traceMatrix( const gsl_matrix* src ){
  int size = GSL_MIN( static_cast<int>(src->size1), static_cast<int>(src->size2) );
  double sum = 0.0;
  for( int i = 0; i < size; ++i ){
    sum += gsl_matrix_get( src, i, i );
  }
  return sum;
}

void normalizeMatrix( gsl_matrix* dest, const gsl_matrix* src ){
  assert( dest->size1 == src->size1 );
  assert( dest->size2 == src->size2 );
  double norm = normMatrix( src );
  gsl_matrix_memcpy( dest, src );
  if( norm > 0.0 ){
    gsl_matrix_scale( dest, 1.0 / norm );
  }else{
    // cerr << __PRETTY_FUNCTION__ << ": norm is zero: " << norm << endl;
    // gsl_matrix_fprintf( stderr, src, "%5.8f " );
    gsl_matrix_set_zero( dest );
  }
  return;
}

double innerProduct( const gsl_matrix* m1, const gsl_matrix* m2 ){
  assert( m1->size1 == m2->size1 );
  assert( m1->size2 == m2->size2 );
  double result = 0.0;
  for( int i = 0; i < static_cast<int>(m1->size1); ++i ){
    for( int j = 0; j < static_cast<int>(m1->size2); ++j ){
      result += gsl_matrix_get( m1, i, j ) * gsl_matrix_get( m2, i, j );
    }
  }
  return result;
}

void invertMatrix( gsl_matrix* dest, const gsl_matrix* src ){
  int signum;
  gsl_matrix* tmp = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_permutation* perm = gsl_permutation_alloc( src->size1 );
  gsl_matrix_memcpy( tmp, src );
  gsl_linalg_LU_decomp( tmp, perm, &signum );
  gsl_linalg_LU_invert( tmp, perm, dest );
  gsl_permutation_free( perm );
  gsl_matrix_free( tmp );
  return;
}

void sqrtMatrix( gsl_matrix* dest, const gsl_matrix* src ){
  gsl_matrix* prevX = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* prevY = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* invX = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* invY = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* XX = gsl_matrix_alloc( src->size1, src->size2 );

  double threshold;
  const double matrixSize = (double)(src->size1 * src->size2 );
  int loop = 0;
  // initialization
  gsl_matrix_memcpy( prevX, src );
  gsl_matrix_set_identity( prevY );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, prevX, prevX, 0.0, XX );
  gsl_matrix_sub( XX, prevX );
  threshold = normMatrix( XX );
  while( (threshold > _EPS) && (loop < _MAXITERATE) ){
    if((detMatrix(prevX) > 0.0) && (detMatrix(prevY) > 0.0)){
      invertMatrix( invX, prevX );
      invertMatrix( invY, prevY );
      //    gsl_matrix_memcpy(tmp, prevX);
      gsl_matrix_add( prevX, invY );
      gsl_matrix_scale( prevX, 0.5 );
      gsl_matrix_add( prevY, invX );
      gsl_matrix_scale( prevY, 0.5 );
      
      // update threshold
      gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, prevX, prevX, 0.0, XX );
      gsl_matrix_sub( XX, src );
      // threshold = normMatrix( XX ) / matrixSize;
      threshold = normMatrix(XX);                    //Bi CK
      ++loop;
    }//if
    else 
      threshold = -1.0; //break
  }//while
  // store result
  gsl_matrix_memcpy( dest, prevX );

  gsl_matrix_free( prevX );
  gsl_matrix_free( prevY );
  gsl_matrix_free( invX );
  gsl_matrix_free( invY );
  gsl_matrix_free( XX );

  return;
}

void logMatrix( gsl_matrix* dest, const gsl_matrix* src ){
  double threshold;
  int k,i;
  int iterate;

  gsl_matrix* I = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* prevA = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* nextA = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* TMP = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* prevZ = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* nextZ = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* X = gsl_matrix_alloc( src->size1, src->size2 );
  
  // initialize 
  gsl_matrix_set_identity( I );
  gsl_matrix_memcpy( prevA, src );
  gsl_matrix_memcpy( TMP, src );
  gsl_matrix_sub( TMP, I );
  k = 0;
  threshold = normMatrix( TMP );

  iterate = 0;
  while( threshold > 0.5 && iterate < _MAXITERATE ){
    sqrtMatrix( nextA, prevA );
    k++;
    gsl_matrix_swap( prevA, nextA );
    gsl_matrix_memcpy( TMP, prevA );
    gsl_matrix_sub( TMP, I );
    //    if (threshold < normMatrix(TMP))
    //      threshold = 0.4;
    //    else
    threshold = normMatrix( TMP );
    ++iterate;
  }
  gsl_matrix_set_identity( nextA );
  gsl_matrix_sub( nextA, prevA );
  gsl_matrix_swap( nextA, prevA );

  gsl_matrix_memcpy( prevZ, prevA );
  gsl_matrix_memcpy( X, prevA );
  i = 1;
  threshold = normMatrix( prevZ );
  iterate = 0;
  while( threshold > _EPS && iterate < _MAXITERATE ){
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, prevZ, prevA, 0.0, nextZ );
    gsl_matrix_swap( nextZ, prevZ );
    gsl_matrix_memcpy( TMP, prevZ );
    i++;
    gsl_matrix_scale( TMP, 1.0 / i );
    gsl_matrix_add( X, TMP );
    //    if(threshold < normMatrix(prevZ))
    //      threshold = 0.0;
    //    else 
    threshold = normMatrix( prevZ );
    ++iterate;
  }
  gsl_matrix_scale( X, pow( 2, k ));
  
  gsl_matrix_memcpy( dest, X );
  
  // free memory
  gsl_matrix_free( X );
  gsl_matrix_free( prevZ );
  gsl_matrix_free( nextZ );
  gsl_matrix_free( prevA );
  gsl_matrix_free( nextA );
  gsl_matrix_free( TMP );
  gsl_matrix_free( I );
  return;
}

void expMatrix( gsl_matrix* dest, const gsl_matrix* src ){
  double c;
  int i, j,k;
  const int q = 6; // this value is come from the thesis
  gsl_matrix* A = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* D = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* invD = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* N = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* prevX = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* nextX = gsl_matrix_alloc( src->size1, src->size2 );
  //  gsl_matrix* XX = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* TMP = gsl_matrix_alloc( src->size1, src->size2 );
  
  // initialize
  gsl_matrix_memcpy( A, src );
    j = 1 + floor( log( normMatrix( A )));
  if( j < 0.0 ) j = 0.0;
  gsl_matrix_scale( A, pow( 2.0, -j ));
  gsl_matrix_set_identity( D );
  gsl_matrix_set_identity( N );
  gsl_matrix_set_identity( prevX );
  c = 1.0;
  
  for( k = 1; k <= q; ++k ){
    c = c * (q - (double)k + 1.0) / ((double)k * (2.0 * q - (double)k + 1.0));
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, A, prevX, 0.0, nextX );
    gsl_matrix_memcpy( prevX, nextX );
    gsl_matrix_scale( nextX, c );
    gsl_matrix_add( N, nextX );
    gsl_matrix_scale( nextX, pow( -1.0, k));
    gsl_matrix_add( D, nextX );
  }
  
  invertMatrix( invD, D );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, invD, N, 0.0, prevX );
  //  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, prevX, prevX, 0.0, XX );
  gsl_matrix_set_identity( dest );

  for( i = 0; i < (int)pow(2, j); ++i ){
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, dest, prevX, 0.0, TMP );
    gsl_matrix_memcpy( dest,  TMP );
  }
  
  gsl_matrix_free( A );
  gsl_matrix_free( D );
  gsl_matrix_free( invD );
  gsl_matrix_free( N );
  gsl_matrix_free( prevX );
  gsl_matrix_free( nextX );
  // gsl_matrix_free( XX );
  gsl_matrix_free( TMP );
  return;
}

void interpolateMatrix( gsl_matrix* dest, const gsl_matrix* src, double ratio ){
  assert( src->size1 == src->size2 );
  assert( dest->size1 == src->size1 );
  assert( dest->size2 == src->size2 );
  if( ratio < 0.0 ) ratio = 0.0;
  if( ratio > 1.0 ) ratio = 1.0;
  
  int i;
  unsigned int bitFlag;

  // make matrix array
  // array[0] = M^(1/2)
  // array[1] = M^(1/4)
  // ...
  // array[n-1] = M^(1/pow(2, n ))
  // const unsigned int intsize = sizeof( int );
  gsl_matrix* tmp = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix** roots = new gsl_matrix*[_INTSIZE];
  for( i = 0; i < _INTSIZE; ++i ){
    roots[i] = gsl_matrix_alloc( src->size1, src->size2 );
  }
  sqrtMatrix( roots[0], src );
  for( i = 0; i + 1 < _INTSIZE; ++i ){
    sqrtMatrix( roots[i+1], roots[i] );
  }

  bitFlag = static_cast<unsigned int>(floor( ratio * static_cast<double>(0xffffffff)));
  gsl_matrix_set_identity( dest );
  for( i = _INTSIZE - 1; i >= 0; i-- ){
    if( (bitFlag & 0x00000001) != 0 ){
      gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, dest, roots[i], 0.0, tmp );
      gsl_matrix_memcpy( dest, tmp);
    }
    bitFlag = bitFlag >> 1;
  }
  for( i = 0; i < _INTSIZE; ++i ){
    gsl_matrix_free( roots[i] );
  }
  delete roots;
  
  return;
}


void logMatrix_sym( gsl_matrix* dest, const gsl_matrix* src ){
  const char funcname[] = "logMatrix_sym: ";
  if( src->size1 != src->size2 ){
    fprintf( stderr, "%s src is not sym matrix\n", funcname );
    return;
  } 
  int i;
  double v;
  gsl_matrix* eigVec = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_vector* eigVal = gsl_vector_alloc( src->size1 );
  eigenSolve_sym( eigVal, eigVec, src, false );
  for( i = 0; i < (int)eigVal->size; ++i ){
    v = gsl_vector_get( eigVal, i );
    if( v > 0.0 ){
      gsl_vector_set( eigVal, i, log(v) );
    }else{
      //cerr << __PRETTY_FUNCTION__ << ": argument is zero" << endl;
      gsl_vector_set( eigVal, i, 0.0 );
    }
  }
  invEigenSolve_sym( dest, eigVal, eigVec );
  gsl_matrix_free( eigVec );
  gsl_vector_free( eigVal );
  return;
}

void expMatrix_sym( gsl_matrix* dest, const gsl_matrix* src ){
  const char funcname[] = "expMatrix_sym: ";
  if( src->size1 != src->size2 ){
    fprintf( stderr, "%s src is not sym matrix\n", funcname );
    return;
  }
  int i;
  double v;
  gsl_matrix* eigVec = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_vector* eigVal = gsl_vector_alloc( src->size1 );
  eigenSolve_sym( eigVal, eigVec, src, false );
  for( i = 0; i < (int)eigVal->size; ++i ){
    v = gsl_vector_get( eigVal, i );
    if( v != -DBL_MAX ){
      gsl_vector_set( eigVal, i, exp(v) );
    }else{
      gsl_vector_set( eigVal, i, DBL_MIN );
    }
  }
  invEigenSolve_sym( dest, eigVal, eigVec );
  gsl_matrix_free( eigVec );
  gsl_vector_free( eigVal );
  return;
}

void isotropicMatrix_sym( gsl_matrix* dest, const gsl_matrix* src ){
  assert( dest->size1 == dest->size2 );
  assert( src->size1 == src->size2 );
  assert( src->size1 == dest->size1 );
  
  double trace = traceMatrix( src );
  gsl_matrix_set_identity( dest );
  gsl_matrix_scale( dest, trace / static_cast<double>(dest->size1) );
  
  return;
}

void deviatoricMatrix_sym( gsl_matrix* dest, const gsl_matrix* src ){
  assert( dest->size1 == dest->size2 );
  assert( src->size1 == src->size2 );
  assert( dest->size1 == src->size1 );
  gsl_matrix* isotropic = gsl_matrix_alloc( src->size1, src->size2 );
  isotropicMatrix_sym( isotropic, src );
  gsl_matrix_memcpy( dest, src );
  gsl_matrix_sub( dest, isotropic );
  gsl_matrix_free( isotropic );
  return;
}

void eigenSolve_sym( gsl_vector* eigVal, gsl_matrix* eigVec, const gsl_matrix* src, const bool zeroCheck ){
//calculate the eigenvalue(eigVal) and eigenvector(eigVec) of matrix(src)

  const char funcname[] = "eigenSolve_sym: ";

  //check whether the matrix is a square matrix
  if( src->size1 != src->size2 ){
    fprintf( stderr, "ERROR: %s it is not symmetric matrix\n", funcname );
    return;
  }

  int i, j, signum;
  gsl_matrix* matrix = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_matrix* vec = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_vector* val = gsl_vector_alloc( src->size1 );
  gsl_permutation* perm = gsl_permutation_alloc( src->size1 );
  gsl_eigen_symmv_workspace* workSpace = gsl_eigen_symmv_alloc( src->size1 );;
  gsl_matrix_memcpy( matrix, src );                                           //copy matrix from "src" to "matrix"
  gsl_eigen_symmv( matrix, val, vec, workSpace );                             //calculate the eigenvalue and eigenvectro of "matrix",
                                                                              //they will be saved in "val" and "vec",respectively.
                                                                              //where, "workSpace" is a temporary space for computation

  gsl_sort_vector_index( perm, val );                                         //sort the eigenvalue by asceding
  gsl_permutation_reverse( perm );                                            //sort the eigenvalue by descending
  
  if( zeroCheck && (gsl_vector_get( val, gsl_permutation_get( perm, (int)perm->size - 1 ) ) <= 0.0 )){
    gsl_matrix_set_identity( vec );
    gsl_vector_set_all( val, 0.0 );
    gsl_permutation_init( perm );
  }
  
  /*
  for( i = 0; i < (int)val->size; ++i ){
    if( zeroCheck && (gsl_vector_get( val, i ) <= 0.0)) gsl_vector_set( val, i, _EPS );
  }
  */
  gsl_sort_vector_index( perm, val );                                         //sort the eigenvalue by asceding
  gsl_permutation_reverse( perm );                                            //sort the eigenvalue by descendin
  //for( i = (int)src->size1 - 1; i >= 0; --i ){
  for( i = 0; i < (int)src->size1; ++i ){                                     //copy the eigenvalue from "val" to "eigVal" by sorting
    gsl_vector_set( eigVal, i, gsl_vector_get( val, gsl_permutation_get( perm, i )) );
    for( j = 0; j < (int)src->size2; ++j ){                                   //copy the eigenvector from "vec" to "eigVec" by sorting
      gsl_matrix_set( eigVec, j, i, gsl_matrix_get( vec, j, gsl_permutation_get( perm, i )));
    }
  }

  gsl_matrix_memcpy( matrix, eigVec );                                        //copy matrix from "eigVec" to "matrix"
  gsl_linalg_LU_decomp( matrix, perm, &signum );                              //decompose "matrix", and save the result into "matrix" and "perm"
  if( gsl_linalg_LU_sgndet( matrix, signum ) < 0 ){                           //calculate the sign of the "matrix", and save it into "signum"
    for( i = 0; i < (int)eigVec->size2; ++i ){
      gsl_matrix_set( eigVec, i, 2, gsl_matrix_get( eigVec, i, 2 ) * -1.0 );
    }
  }

  gsl_matrix_free( matrix );
  gsl_matrix_free( vec );
  gsl_vector_free( val );
  gsl_permutation_free( perm );
  gsl_eigen_symmv_free( workSpace );
  return;
}

void invEigenSolve_sym( gsl_matrix* dest, const gsl_vector* eigVal, const gsl_matrix* eigVec ){
  //calculate the matrix of tensor
  int i;
  gsl_matrix* temp = gsl_matrix_alloc( eigVec->size1, eigVec->size2 );
  gsl_matrix* diag = gsl_matrix_alloc( eigVec->size1, eigVec->size2 );
  gsl_matrix_set_zero( diag );                             
  for( i = 0; i < (int)eigVal->size; ++i ){
    gsl_matrix_set( diag, i, i, gsl_vector_get( eigVal, i ));
  }

  //calculate the matrix of tensor by using (eigenvector)*(eigenvalue)*(inverse of (eigenvector))
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, eigVec, diag, 0.0, temp );
  gsl_matrix_memcpy( diag, temp );

    gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, diag, eigVec, 0.0, temp ); //CblasTrans means inverse of eigVec
  //calculate the inverse matrix of eigVec
  //  int size = eigVal->size;
  // gsl_permutation*perm = gsl_permutation_alloc(size);
  // gsl_matrix*inverse = gsl_matrix_alloc(eigVec->size1, eigVec->size2);
  // gsl_matrix* temp_eigVal = gsl_matrix_alloc( eigVec->size1, eigVec->size2 );
  // gsl_matrix_memcpy(temp_eigVal, eigVec);
  // int s=0;
  // gsl_linalg_LU_decomp(temp_eigVal, perm, &s);
  // gsl_linalg_LU_invert(temp_eigVal, perm, inverse);
  // calculate the inverse matrix of eigVec

  // gsl_matrix* temp1 = gsl_matrix_alloc( eigVec->size1, eigVec->size2 );
  // gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, diag, inverse, 0.0, temp1 );
  //calculate the matrix of tensor by using (eigenvector)*(eigenvalue)*(inverse of (eigenvector))

  gsl_matrix_memcpy( dest, temp );
  gsl_matrix_free( temp );
  // gsl_matrix_free( temp1 );
  gsl_matrix_free( diag );
  // gsl_matrix_free( temp_eigVal );
  // gsl_matrix_free( inverse );
  // gsl_permutation_free(perm);

  return;
}

double detMatrix( const gsl_matrix* src ){
  gsl_matrix* matrix = gsl_matrix_alloc( src->size1, src->size2 );
  gsl_permutation* perm = gsl_permutation_alloc( src->size1 );
  int signum;
  gsl_matrix_memcpy( matrix, src );
  gsl_linalg_LU_decomp( matrix, perm, &signum );
  double det = gsl_linalg_LU_det( matrix, signum );
  gsl_matrix_free( matrix );
  gsl_permutation_free( perm );
  return det;
}

void rotationMatrix( gsl_matrix* dest, const gsl_vector* axis, const double degree ){
  if( fabs(degree) < _EPS ){
    gsl_matrix_set_identity( dest );
    return;
  }
  gsl_vector* vec = gsl_vector_alloc( axis->size );
  gsl_vector_memcpy( vec, axis );
  gsl_vector_scale( vec, degree / gsl_blas_dnrm2( vec ));
  double EM3[3];
  double R[4][4];
  int i, j;
  for( i = 0; i < 3; ++i ){
    EM3[i] = gsl_vector_get( vec, i );
  }
  EM3_To_R( EM3, R );
  for( i = 0; i < 3; ++i ){
    for( j = 0; j < 3; ++j ){
      gsl_matrix_set( dest, i, j, R[i][j] );
    }
  }
  gsl_vector_free( vec );
  return;
}

void getRotMatrix( gsl_matrix* dest, const gsl_matrix* src1, const gsl_matrix* src2 ){
  gsl_matrix* invSrc1 = gsl_matrix_alloc( src1->size1, src1->size2 );
  invertMatrix( invSrc1, src1 );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, src2, invSrc1, 0.0, dest );
  gsl_matrix_free( invSrc1 );
  return;
}

void interpolateRotMatrix( gsl_matrix* dest, const gsl_matrix* src, const double param ){
  interpolateMatrix( dest, src, param );
  return;
}


void expMap( gsl_matrix* rotMatrix, const gsl_vector* axis, const double radian ){
  // check arguments
  if( rotMatrix->size1 != 3 || rotMatrix->size2 != 3 ){
    //fprintf( stderr, "%s: invalid size ", __PRETTY_FUNCTION__ );
    return;
  }
  if( axis->size != 3 ){
    //fprintf( stderr, "%s: invalid size ", __PRETTY_FUNCTION__ );
    return;
  }
  
  // if angle = 0 -> rotMatrix = identity
  if( radian == 0.0 ){
    gsl_matrix_set_identity( rotMatrix );
    return;
  }
  
  // variables
  double AXIS[3];
  double ROTMATRIX[4][4];
  
  for( int i = 0; i < 3; ++i ){
    AXIS[i] = gsl_vector_get( axis, i );
    AXIS[i] *= radian / gsl_blas_dnrm2( axis );
  }
  
  EM3_To_R( AXIS, ROTMATRIX );
  
  for( int i = 0; i < 3; ++i ){
    for( int j = 0; j < 3; ++j ){
      gsl_matrix_set( rotMatrix, i, j, ROTMATRIX[i][j] );
    }
  }
  
  return;
}
