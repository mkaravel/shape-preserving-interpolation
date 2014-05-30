// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_MATH_IPOWER_HPP
#define SHAPE_PRESERVING_UTIL_MATH_IPOWER_HPP


namespace shape_preserving
{

namespace util { namespace math
{


template <typename NT>
inline NT ipower(const NT& x, unsigned int n)
{
    NT y = x;
    NT z(1);

    //  z = x;
    unsigned int m = n;
    while (n > 1)
    {
        m = m / 2;
        if ( n > 2 * m )
        {
            z = z * y;
        }
        y = y * y;
        n = m;
    }

    z = z * y;
    return z;
}


}} // namespace util::math

} // namespace shape_preserving


#endif // SHAPE_PRESERVING_UTIL_MATH_IPOWER_HPP
