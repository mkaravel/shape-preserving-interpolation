// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_MATH_BERNSTEIN_3_HPP
#define SHAPE_PRESERVING_UTIL_MATH_BERNSTEIN_3_HPP


namespace shape_preserving
{

namespace util { namespace math
{

template <class NT_t>
struct bernstein_3
{
    typedef NT_t  NT;

    int degree() const { return 3; }

    NT value(int i, const NT& t) const
    {
        switch (i)
        {
        case 0:
            {
                NT tt = NT(1) - t;
                return tt * tt * tt;
                //	break;
            }
        case 1:
            {
                NT tt = NT(1) - t;
                return NT(3) * tt * tt * t;
                //      break;
            }
        case 2:
            return NT(3) * (NT(1) - t) * t * t;
            //      value = 3 * (1 - t) * t * t;
            //      break;
        case 3:
            return t * t * t;
            //      value = t * t * t;
            //      break;
        default:
            return 0;
            //      value = 0;
            //      break;
        }

        return 0;
  }

  NT der1(int i, const NT& t) const
  {
      switch (i)
      {
      case 0:
          {
              NT tt = NT(1) - t;
              return NT(-3) * tt * tt;
              //    break;
          }
      case 1:
          return NT(3) * (NT(1) - t) * (NT(1) - NT(3) * t);
          //      break;
      case 2:
          return NT(3) * t * (NT(2) - NT(3) * t);
          //      break;
      case 3:
          return NT(3) * t * t;
          //      break;
      default:
          return 0;
          //      break;
      }
      return 0;
  }

  NT der2(int i, const NT& t) const
  {
      switch (i)
      {
      case 0:
          return NT(6) * (NT(1) - t);
          //      break;
      case 1:
          return NT(-6) * (NT(2) - NT(3) * t);
          //      break;
      case 2:
          return NT(6) * (NT(1) - NT(3) * t);
          //      break;
      case 3:
          return NT(6) * t;
          //      break;
      default:
          return 0;
          break;
      }

      return 0;
  }

  NT der3(int i, const NT&) const
  {
      switch (i)
      {
      case 0:
          return -6;
          //      break;
      case 1:
          return 18;
          //      break;
      case 2:
          return -18;
          //      break;
      case 3:
          return 6;
          //      break;
      default:
          return 0;
          //      break;
      }

      return 0;
  }
};


}} // namespace util::math

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_MATH_BERNSTEIN_3_HPP
