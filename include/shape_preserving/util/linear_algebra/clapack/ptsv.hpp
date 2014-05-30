#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTSV_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTSV_HPP

#include <algorithm>

#include <shape_preserving/util/linear_algebra/clapack/xerbla.hpp>
#include <shape_preserving/util/linear_algebra/clapack/pttrf.hpp>
#include <shape_preserving/util/linear_algebra/clapack/pttrs.hpp>

LINEAR_ALGEBRA_BEGIN_NAMESPACE

template <typename NT>
inline int ptsv(const int& n, const int& nrhs,  NT *d, 
                NT *e, NT *b, const int& ldb, int& info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTSV computes the solution to a real system of linear equations   
    A*X = B, where A is an N-by-N symmetric positive definite tridiagonal 
  
    matrix, and X and B are N-by-NRHS matrices.   

    A is factored as A = L*D*L**T, and the factored form of A is then   
    used to solve the system of equations.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.  On exit, the n diagonal elements of the diagonal matrix   
            D from the factorization A = L*D*L**T.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A.  On exit, the (n-1) subdiagonal elements of the   
            unit bidiagonal factor L from the L*D*L**T factorization of   
            A.  (E can also be regarded as the superdiagonal of the unit 
  
            bidiagonal factor U from the U**T*D*U factorization of A.)   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the solution has not been   
                  computed.  The factorization has not been completed   
                  unless i = N.   

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
    xerbla("PTSV", -info);
    return 0;
  }

  /*     Compute the L*D*L' (or U'*D*U) factorization of A. */

  pttrf(n, &d[0], &e[0], info);
  if (info == 0) {
    /*        Solve the system A*X = B, overwriting B with X. */

    pttrs(n, nrhs, &d[0], &e[0], &b[0], ldb, info);
  }
  return 0;

/*     End of PTSV */

} /* ptsv */

LINEAR_ALGEBRA_END_NAMESPACE

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_PTSV_HPP
