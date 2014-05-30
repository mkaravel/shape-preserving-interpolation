// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_GCTSV_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_GCTSV_HPP

#include <cassert>

#include <shape_preserving/util/linear_algebra/clapack/gtsv.hpp>
#include <shape_preserving/util/linear_algebra/ctsvsm.hpp>


namespace shape_preserving
{

namespace util { namespace linear_algebra
{


// Solves the equation
//
//    A * X = B
//
// where A is an n by n cyclic-tridiagonal matrix, by writting
// A as T+u.v^t, where T is tridiagonal, and using the
// Sherman-Morrison formula
//
// The input is an in gtsv; the top-right element of A is stored at
// du[n-1], whereas the bottom-left element of A is stored at dl[n-1]

template <typename NT>
inline int gctsv(const int& n, const int& nrhs, NT *dl,
                 NT *d, NT* du, NT *b, const int& ldb, int& info)
{
  // dl and du are overwritten. So even if the upper and lower
  // diagonal elements are equal, we must pass different copies of
  // them to the subroutine
  assert( du != dl );

  // Update diagonal elements
  d[0] -= du[n - 1];
  d[n - 1] -= dl[n - 1];

  // Allocate and evaluate local_b
  int newnrhs = nrhs + 1;
  NT *local_b = new NT[newnrhs * n];
  sm_rhs(n, nrhs, b, du[n - 1], dl[n - 1], local_b);

  // Solve the linear systems
  (void) clapack::gtsv(n, newnrhs, dl, d, du, local_b, n, info);

  if (info == 0) {
    sm_solution(n, nrhs, local_b, b);
  }

  delete []local_b;

  return info;
}

}} // namespace util::linear_algebra

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_GCTSV_HPP
