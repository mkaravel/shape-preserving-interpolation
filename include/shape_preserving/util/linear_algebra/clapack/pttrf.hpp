#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTTRF_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTTRF_HPP

#include <shape_preserving/util/linear_algebra/clapack/xerbla.hpp>

namespace shape_preserving
{

namespace util { namespace linear_algebra
{

namespace clapack
{

template <typename NT>
inline int pttrf(const int& n, NT *d, NT *e, int& info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTTRF computes the factorization of a real symmetric positive   
    definite tridiagonal matrix A.   

    If the subdiagonal elements of A are supplied in the array E, the   
    factorization has the form A = L*D*L**T, where D is diagonal and L   
    is unit lower bidiagonal; if the superdiagonal elements of A are   
    supplied, it has the form A = U**T*D*U, where U is unit upper   
    bidiagonal.  (The two forms are equivalent if A is real.)   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.  On exit, the n diagonal elements of the diagonal matrix   
            D from the L*D*L**T factorization of A.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) off-diagonal elements of the tridiagonal 
  
            matrix A.   
            On exit, the (n-1) off-diagonal elements of the unit   
            bidiagonal factor L or U from the factorization of A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite; if i < N, the factorization could   
                  not be completed, while if i = N, the factorization was 
  
                  completed, but D(N) = 0.   

    ===================================================================== 

    November 2007:
    Adapted by Menelaos Karavelas to be a template function

    =====================================================================  

*/
  static int i;
  static NT di, ei;

  info = 0;
  if (n < 0) {
    info = -1;
    xerbla("PTTRF", -info);
    return 0;
  }

  /*     Quick return if possible */

  if (n == 0) {	return 0; }

  /*     Compute the L*D*L' (or U'*D*U) factorization of A. */

  for (i = 1; i <= n - 1; ++i) {

    /*      Drop out of the loop if d(i) <= 0: the matrix is not
	    positive definite. */

    di = d[i-1];
    if (di <= 0) {
      info = i;
      return 0;
    }

    /*        Solve for e(i) and d(i+1). */

    ei = e[i-1];
    e[i-1] = ei / di;
    d[i] -= e[i-1] * ei;
  }

  /*     Check d(n) for positive definiteness. */

  i = n;
  if (d[i-1] > 0) { return 0; }

  info = i;

  return 0;

/*     End of PTTRF */

} /* pttrf */

} // namespace clapack

}} // namespace util::linear_algebra

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTTRF_HPP
