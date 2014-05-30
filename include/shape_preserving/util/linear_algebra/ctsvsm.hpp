// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CTSVSM_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CTSVSM_HPP

namespace shape_preserving
{

namespace util { namespace linear_algebra
{


template <typename NT>
inline void sm_rhs(const int& n, const int& nrhs, const NT* bc,
                   const NT& ur, const NT& bl, NT* b)
{
    // evaluate new RHS
    int nrhs_x_n = nrhs * n;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < nrhs; ++j)
        {
            b[i + j * n] = bc[i + j * n];
        }

        if (i == 0)
        {
            b[nrhs_x_n] = ur;
        }
        else if (i == n - 1)
        {
            b[nrhs_x_n + n - 1] = bl;
        }
        else
        {
            b[i + nrhs_x_n] = NT(0);
        }
    }
}


template <typename NT>
inline void sm_solution(const int& n, const int& nrhs, const NT* b, NT* bc)
{
    // compute solution for original cyclic-tridiagonal system
    int nrhs_x_n = nrhs * n;

    NT alpha = NT(-1) / (NT(1) + b[nrhs_x_n] + b[nrhs_x_n + n - 1]);

    for (int j = 0; j < nrhs; ++j)
    {
        int jn = j * n;
        NT alpha_dot = alpha * (b[jn] + b[jn + n - 1]);

        for (int i = 0; i < n; ++i)
        {
            bc[i + jn] = b[i + jn] + alpha_dot * b[i + nrhs_x_n];
        }
    }
}

}} // namespace util::linear_algebra

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CTSVSM_HPP
