// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_MATH_POLYNOMIAL_VECTOR_3_HPP
#define SHAPE_PRESERVING_UTIL_MATH_POLYNOMIAL_VECTOR_3_HPP

#include <iostream>
#include <shape_preserving/geometries/vector_3.h>


namespace shape_preserving
{

namespace util { namespace math
{


template<class P>
class polynomial_vector_3
{
public:
  typedef P polynomial;

  typedef polynomial_vector_3<polynomial>  Self;
  typedef typename polynomial::NT          NT;
  typedef polynomial                       Coordinate;

private:
  polynomial coord[3];

public:
  polynomial_vector_3() {};
  polynomial_vector_3(const polynomial& v1,
		      const polynomial& v2,
		      const polynomial& v3)
  { 
    coord[0] = v1; coord[1] = v2; coord[2] = v3;
  }
  
  inline const polynomial& operator[](unsigned int i) const {
    return coord[i];
  }

  inline polynomial& operator[](unsigned int i) {
    return coord[i];
  }

  bool operator==(const Self& q) const
  {
    return (coord[0] == q[0] && coord[1] == q[1] && coord[2] == q[2]);
  }

  inline Self& operator+() { return (*this); }

  inline Self operator-() const {
    return Self(-coord[0], -coord[1], -coord[2]);
  }

  Self operator-(const Self& q) const
  {
    Self result(*this);
    result -= q;
    return result;
  }

  Self operator+(const Self& q) const
  {
    Self result(*this);
    result += q;
    return result;
  }

  Self& operator-=(const Self& q)
  {
    coord[0] -= q[0];
    coord[1] -= q[1];
    coord[2] -= q[2];

    return *this;
  }

  Self& operator+=(const Self& q)
  {
    coord[0] += q[0];
    coord[1] += q[1];
    coord[2] += q[2];

    return *this;
  }

  Self operator*(const NT& a) const
  {
    Self result(*this);
    result *= a;
    return result;
  }

  Self& operator*=(const NT& a)
  {
    coord[0] *= a;
    coord[1] *= a;
    coord[2] *= a; 

    return *this;
  }

  polynomial operator*(const Self& q) const
  {
    return (coord[0] * q[0] + coord[1] * q[1] + coord[2] * q[2]);
  }

  Self operator/(const NT& a) const
  {
    Self result(*this);
    result /= a;
    return result;
  }

  Self& operator/=(const NT& a)
  {
    coord[0] /= a;
    coord[1] /= a;
    coord[2] /= a; 

    return *this;
  }

  Self cross(const Self& q) const
  {
    Self result;

    result.coord[0] = coord[1] * q[2] - coord[2] * q[1];
    result.coord[1] = coord[2] * q[0] - coord[0] * q[2];
    result.coord[2] = coord[0] * q[1] - coord[1] * q[0];

    return result;
  }

  inline void degree(int d)
  {
    for (int i = 0; i < 3; i++) { coord[i].degree(d); }
  }

  inline Self grad() const
  {
    return Self(coord[0].derivative(),
		coord[1].derivative(),
		coord[2].derivative());
  }


  Vector_3<NT> value(const NT& x)
  {
    Vector_3<NT> val;

    for (int i = 0; i < 3; i++) {
      val[i] = coord[i].value(x);
    }

    return val;
  }

};

template<typename NT, class P>
polynomial_vector_3<P>
cross_product(const polynomial_vector_3<P>& pv, const Vector_3<NT>& v)
{
  polynomial_vector_3<P> result;

  result[0] = pv[1] * v[2] - pv[2] * v[1];
  result[1] = pv[2] * v[0] - pv[0] * v[2];
  result[2] = pv[0] * v[1] - pv[1] * v[0];

  return result;
}


template<class P>
polynomial_vector_3<P>
operator*(const typename P::NT& a, const polynomial_vector_3<P>& v)
{
  polynomial_vector_3<P> result = v;
  result *= a;
  return result;
}


template<class P>
std::ostream& operator<<(std::ostream &o, const polynomial_vector_3<P>& v)
{
  o << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return o;
}


}} // namespace util::math

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_MATH_POLYNOMIAL_VECTOR_3_HPP
