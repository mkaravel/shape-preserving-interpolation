#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTTRS_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTTRS_HPP

#include <algorithm>

#include <shape_preserving/util/linear_algebra/clapack/xerbla.hpp>

namespace shape_preserving
{

namespace util { namespace linear_algebra
{

namespace clapack
{

template <typename NT>
inline int pttrs(const int& n, const int& nrhs, NT *d, 
                 NT *e, NT *b, const int& ldb, int& info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTTRS solves a system of linear equations A * X = B with a   
    symmetric positive definite tridiagonal matrix A using the   
    factorization A = L*D*L**T or A = U**T*D*U computed by DPTTRF.   
    (The two forms are equivalent if A is real.)   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the tridiagonal matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization computed by DPTTRF.   

    E       (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) off-diagonal elements of the unit bidiagonal factor 
  
            U or L from the factorization computed by DPTTRF.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 

    November 2007:
    Adapted by Menelaos Karavelas to be a template function

    =====================================================================  

*/

  info = 0;
  if (n < 0) {
    info = -1;
  } else if (nrhs < 0) {
    info = -2;
  } else if (ldb < (std::max)(1,n)) {
    info = -6;
  }
  if (info != 0) {
    xerbla("PTTRS", -info);
    return 0;
  }

  /*     Quick return if possible */

  if ( n == 0 ) { return 0; }

  /*     Solve A * X = B using the factorization A = L*D*L',   
	 overwriting each right hand side vector with its solution. */

  for (int j = 1; j <= nrhs; ++j) {
    /*        Solve L * x = b. */
    for (int i = 2; i <= n; ++i) {
      int k = i - 1 + (j - 1) * ldb;
      int k1 = k - 1;
      b[k] -= b[k1] * e[i-2];
    }

    /*        Solve D * L' * x = b. */

    b[n - 1 + (j - 1) * ldb] /= d[n - 1];
    for (int i = n - 1; i >= 1; --i) {
      int k = i + (j - 1) * ldb;
      int k1 = k - 1;
      b[k1] = b[k1] / d[i-1] - b[k] * e[i-1];
    }
  }
  
  return 0;

  /*     End of PTTRS */

} /* pttrs */

} // namespace clapack

}} // namespace util::linear_algebra

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTTRS_HPP
