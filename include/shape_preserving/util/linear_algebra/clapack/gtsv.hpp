#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_GTSV_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_GTSV_HPP

#include <cassert>
#include <cstdlib>
#include <cmath>

#include <algorithm>

#include <shape_preserving/util/linear_algebra/clapack/xerbla.hpp>


namespace shape_preserving
{

namespace util { namespace linear_algebra
{

namespace clapack
{

template <typename NT>
inline int gtsv(const int& n, const int& nrhs, NT *dl,
                NT *d, NT *du, NT *b, const int& ldb, int& info)
{
  // dl and du are overwritten. So even if the upper and lower
  // diagonal elements are equal, we must pass different copies of
  // them to the subroutine
  assert( du != dl );

  /* System generated locals */
  NT d1, d2;


#define b_ref(a_1,a_2) b[(a_2) * ldb + a_1]


/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1999   


    Purpose   
    =======   

    DGTSV  solves the equation   

       A*X = B,   

    where A is an n by n tridiagonal matrix, by Gaussian elimination with   
    partial pivoting.   

    Note that the equation  A'*X = B  may be solved by interchanging the   
    order of the arguments DU and DL.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DL must contain the (n-1) sub-diagonal elements of   
            A.   

            On exit, DL is overwritten by the (n-2) elements of the   
            second super-diagonal of the upper triangular matrix U from   
            the LU factorization of A, in DL(1), ..., DL(n-2).   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, D must contain the diagonal elements of A.   

            On exit, D is overwritten by the n diagonal elements of U.   

    DU      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DU must contain the (n-1) super-diagonal elements   
            of A.   

            On exit, DU is overwritten by the (n-1) elements of the first   
            super-diagonal of U.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N by NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N by NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, U(i,i) is exactly zero, and the solution   
                 has not been computed.  The factorization has not been   
                 completed unless i = N.   

    =====================================================================   

    November 2007:
    Adapted by Menelaos Karavelas to be a template function

    =====================================================================
*/

  // Parameter adjustments
  --dl;
  --d;
  --du;
  b -= 1 + ldb;

  // Function Body
  info = 0;
  if (n < 0) {
    info = -1;
  } else if (nrhs < 0) {
    info = -2;
  } else if (ldb < (std::max)(1,n)) {
    info = -7;
  }
  if (info != 0) {
    xerbla("GTSV", -info);
    return 0;
  }

  if (n == 0) {
    return 0;
  }

  if (nrhs == 1) {
    int i1 = n - 2;
    for (int i = 1; i <= i1; ++i) {
      int ip = i + 1;
      if ((d1 = d[i], std::abs(d1)) >= (d2 = dl[i], std::abs(d2))) {
	// No row interchange required
	if (d[i] != NT(0)) {
	  NT fact = dl[i] / d[i];
	  d[ip] -= fact * du[i];
	  b_ref(ip, 1) = b_ref(ip, 1) - fact * b_ref(i, 1);
	} else {
	  info = i;
	  return 0;
	}
	dl[i] = NT(0);
      } else {
	// Interchange rows I and I+1
	NT fact = d[i] / dl[i];
	d[i] = dl[i];
	NT temp = d[ip];
	d[ip] = du[i] - fact * temp;
	dl[i] = du[ip];
	du[ip] = -fact * dl[i];
	du[i] = temp;
	temp = b_ref(i, 1);
	b_ref(i, 1) = b_ref(ip, 1);
	b_ref(ip, 1) = temp - fact * b_ref(ip, 1);
      }
    }
    if (n > 1) {
      int i = n - 1;
      int ip = n;
      if ((d1 = d[i], std::abs(d1)) >= (d2 = dl[i], std::abs(d2))) {
	if (d[i] != NT(0)) {
	  NT fact = dl[i] / d[i];
	  d[ip] -= fact * du[i];
	  b_ref(ip, 1) = b_ref(ip, 1) - fact * b_ref(i, 1);
	} else {
	  info = i;
	  return 0;
	}
      } else {
	NT fact = d[i] / dl[i];
	d[i] = dl[i];
	NT temp = d[ip];
	d[ip] = du[i] - fact * temp;
	du[i] = temp;
	temp = b_ref(i, 1);
	b_ref(i, 1) = b_ref(ip, 1);
	b_ref(ip, 1) = temp - fact * b_ref(ip, 1);
      }
    }
    if (d[n] == NT(0)) {
      info = n;
      return 0;
    }
  } else {
    int i1 = n - 2;
    for (int i = 1; i <= i1; ++i) {
      int ip = i + 1;
      if ((d1 = d[i], std::abs(d1)) >= (d2 = dl[i], std::abs(d2))) {
	// No row interchange required
	if (d[i] != NT(0)) {
	  NT fact = dl[i] / d[i];
	  d[ip] -= fact * du[i];
	  for (int j = 1; j <= nrhs; ++j) {
	    b_ref(ip, j) = b_ref(ip, j) - fact * b_ref(i, j);
	  }
	} else {
	  info = i;
	  return 0;
	}
	dl[i] = NT(0);
      } else {
	// Interchange rows I and I+1
	NT fact = d[i] / dl[i];
	d[i] = dl[i];
	NT temp = d[ip];
	d[ip] = du[i] - fact * temp;
	dl[i] = du[ip];
	du[ip] = -fact * dl[i];
	du[i] = temp;
	for (int j = 1; j <= nrhs; ++j) {
	  temp = b_ref(i, j);
	  b_ref(i, j) = b_ref(i + 1, j);
	  b_ref(ip, j) = temp - fact * b_ref(ip, j);
	}
      }
    }

    if (n > 1) {
      int i = n - 1;
      int ip = n;
      if ((d1 = d[i], std::abs(d1)) >= (d2 = dl[i], std::abs(d2))) {
	if (d[i] != NT(0)) {
	  NT fact = dl[i] / d[i];
	  d[ip] -= fact * du[i];
	  for (int j = 1; j <= nrhs; ++j) {
	    b_ref(ip, j) = b_ref(ip, j) - fact * b_ref(i, j);
	  }
	} else {
	  info = i;
	  return 0;
	}
      } else {
	NT fact = d[i] / dl[i];
	d[i] = dl[i];
	NT temp = d[ip];
	d[ip] = du[i] - fact * temp;
	du[i] = temp;
	for (int j = 1; j <= nrhs; ++j) {
	  temp = b_ref(i, j);
	  b_ref(i, j) = b_ref(i + 1, j);
	  b_ref(ip, j) = temp - fact * b_ref(ip, j);
	}
      }
    }
    if (d[n] == NT(0)) {
      info = n;
      return 0;
    }
  }

  // Back solve with the matrix U from the factorization
  if (nrhs <= 2) {
    int j = 1;
    while ( true ) {
      b_ref(n, j) = b_ref(n, j) / d[n];
      if (n > 1) {
	int n1 = n - 1;
	b_ref(n1, j) = (b_ref(n1, j) - du[n1] * b_ref(n, j)) / d[n1];
      }
      for (int i = n - 2; i >= 1; --i) {
	b_ref(i, j) = (b_ref(i, j) - du[i] * b_ref(i + 1, j) 
		       - dl[i] * b_ref(i + 2, j)) / d[i];
      }
      if (j >= nrhs) { break; }
      ++j;
    }
  } else {
    for (int j = 1; j <= nrhs; ++j) {
      b_ref(n, j) = b_ref(n, j) / d[n];
      if (n > 1) {
	int n1 = n - 1;
	b_ref(n1, j) = (b_ref(n1, j) - du[n1] * b_ref(n, j)) / d[n1];
      }
      for (int i = n - 2; i >= 1; --i) {
	b_ref(i, j) = (b_ref(i, j) - du[i] * b_ref(i + 1, j) 
		       - dl[i] * b_ref(i + 2, j)) / d[i];
      }
    }
  }

  return 0;

/*     End of GTSV */

} /* gtsv */

#undef b_ref


} // namespace clapack

}} // namespace util::linear_algebra

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_GTSV_HPP
