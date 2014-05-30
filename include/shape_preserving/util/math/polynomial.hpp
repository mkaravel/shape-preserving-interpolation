// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_MATH_POLYNOMIAL_HPP
#define SHAPE_PRESERVING_UTIL_MATH_POLYNOMIAL_HPP

#include <cstddef>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <algorithm>


namespace shape_preserving
{

namespace util { namespace math
{


template <typename NT_t, std::size_t MaxDegree = 10>
class polynomial
{
public:
    typedef NT_t             NT;
    typedef polynomial<NT>   Self;

    static const std::size_t max_degree = MaxDegree;

private:
    static const double SAME_COEF_EPS = 0;
    static const int IOS_PREC = 20;

    int deg;
    NT coefs[MaxDegree];

    bool same_coefs(const Self& q) const
    {
        int i(0);

        while (i <= deg && 
               std::abs(coefs[i] - q.coefs[i]) <= SAME_COEF_EPS *   
               (std::max)(std::abs(coefs[i]), std::abs(q.coefs[i])))
        {
            i++;
        }

        if (i == deg + 1) { return true; }
        return false;
    }

    void recompute_degree()
    {
        // Decreases the degree if last coefs are 0.
        while ( deg > 0 && coefs[deg] == 0 )
        {
            deg--;
        }
    }

    void make_sure_zero(int n)
    {
        // Make sure coefs from degree+1 to n are 0  
        for (int i = degree() + 1 ; i <= n ; ++i)
        {
            coefs[i] = 0;
        }
        if (n > degree()) { deg = n; }
    }

    Self mult(const Self& q) const // without origin
    {
        int n = degree() + q.degree();
        Self result(n);

        for (int i = degree(); i >= 0; --i)
        {
            for (int j = q.degree(); j >= 0; --j)
            {
                result.coefs[i+j] += coefs[i] * q.coef(j);
            }
        }

        result.recompute_degree();
        return result;
    }

public:
    polynomial()
    {
        deg = 0;
        coefs[0] = 0;
    }

    polynomial(int d)
    {
        deg = d;
        for (int i = 0; i <= d; ++i) { coefs[i] = 0; }
    }

    polynomial(const NT& a)
    {
        coefs[0] = a;
        deg = 0;
    }

    polynomial(const NT& a, const NT& b)
    {
        coefs[0] = a;
        coefs[1] = b;
        deg = 1;
    }

    polynomial(const NT a[], int d)
    {
        assert(d < MAXPOLYDEG);
        assert(sizeof(a) / sizeof(a[0]) < d + 1);

        deg = d;
        for (int i = 0; i <= d; i++)
        {
            coefs[i] = a[i];
        }
    }

    polynomial(const Self& p)
    {
        assert(p.deg < MAXPOLYDEG);
    
        deg = p.deg;
        for (int i = 0; i <= deg; i++)
        {
            coefs[i] = p.coefs[i];
        }
    }

    inline int degree() const { return deg; }
    inline void degree(int d) { deg = d;  recompute_degree(); } 

    void coef(int i, const NT& c)
    {
        assert( i < MAXPOLYDEG );
        coefs[i] = c;
        if ( i > deg )
        {
            for(int j = deg+1 ; j < i ; j++)
            {
                coefs[j] = 0;
            }
            deg = i;
        }
        recompute_degree();
    }

    inline const NT& coef(int i) const
    {
        assert( i <= deg );
        return coefs[i];
    }
  
    inline bool is_zero() const { return (deg == 0 && coefs[0] == 0); }

    inline const NT& operator[](int i) const { return coef(i); }
    bool operator==(const Self& q) const
    {
        return (deg == q.deg && same_coefs(q));
    }

    Self operator-()
    {
        Self res;
    
        for (int i = 0; i <= deg; i++)
        {
            res.coefs[i] = -coefs[i];
        }
        res.degree(deg);

        return res;
    }

    Self operator-(const Self& q) const
    {
        Self result = *this;
        result -= q;
        return result;
    }

    Self operator+(const Self& q) const
    {
        Self result = *this;
        result += q;
        return result;
    }

    Self& operator-=(const Self& q)
    {
        make_sure_zero(q.degree());

        for (int i = q.degree(); i >= 0 ; --i)
        {
            coefs[i] -= q.coef(i);
        }

        recompute_degree();
        return *this;
    }

    Self& operator+=(const Self& q)
    {
        make_sure_zero(q.degree());

        for (int i = q.degree(); i >= 0 ; --i)
        {
            coefs[i] += q.coef(i);
        }

        recompute_degree();
        return *this;
    }

    Self operator-(const NT& c) const
    {
        Self res(*this);
        res.coefs[0] -= c;
        return res;
    }

    Self operator+(const NT& c) const
    {
        Self res(*this);
        res.coefs[0] += c;
        return res;
    }

    Self& operator-=(const NT& c)
    {
        coefs[0] -= c;
        return *this;
    }

    Self& operator+=(const NT& c)
    {
        coefs[0] += c;
        return *this;
    }

    Self& operator*=(const NT& a)
    {
        for (int i = 0; i <= deg; ++i)  { coefs[i] *= a; }
        return *this;
    }

    Self operator*(const NT& a) const
    {
        Self result(*this);
        result *= a;
        return result;
    }

    /* polynomial multiplication */
    Self operator*(const Self& q) const
    {
        int n = degree() + q.degree();
        Self result(n);
        int i, j;

        for (i = degree(); i >= 0; --i)
        {
            for (j = q.degree(); j >= 0; --j)
            {
                result.coefs[i+j] += coefs[i] * q.coef(j);
            }
        }

        result.recompute_degree();
        return result;
    }

    Self& operator*=(const Self& q)
    {
        *this = (*this * q);
        return *this;
    }

    /* polynomial division */
    Self operator/(const Self& q) const
    {
        const int n = q.degree();
        Self u = *this;
        int divdeg = degree() - n;
        assert( divdeg >= 0 );
        Self div(divdeg);

        NT v_ninv = NT(1) / q.coefs[n];

        for (int k = divdeg; k >= 0; k--)
        {
            div.coefs[k] = u.coefs[n + k] * v_ninv;
            for (int j = n + k - 1; j >= k; j--)
            {
                u.coefs[j] -= div.coefs[k] * q.coefs[j - k];
            }
        }

        return (div);
    }
    
    Self& operator/=(const Self& q)
    {
        *this = (*this / q);
        return *this;
    }

    /* polynomial modulo */
    Self operator%(const Self& q) const
    {
        const int n = q.degree();
        int divdeg = degree() - n;
        assert( divdeg >= 0 );
        Self modulo;

        if (n == 0) return (modulo);

        modulo = *this;
        Self div(divdeg);

        NT v_ninv = NT(1) / q.coefs[n];
  
        for (int k = divdeg; k >= 0; k--)
        {
            div.coefs[k] = modulo.coefs[n + k] * v_ninv;
            for (int j = n + k - 1; j >= k; j--)
            {
                modulo.coefs[j] -= div.coefs[k] * q.coefs[j - k];
            }
        }

        modulo.deg = n - 1;
        return modulo;
    }

    Self& operator%=(const Self& q)
    {
        *this = (*this % q);
        return *this;
    }

    /* Derivative polynomial */
    Self derivative() const
    {
        int n = degree() - 1;
        Self result(n);
  
        while (n >= 0)
        {
            result.coefs[n] = NT(n+1) * coefs[n+1];
            n--;
        }

        result.recompute_degree();
        return result;
    }

    /* integral of the polynomial from a to b */
    NT integral(const NT& a, const NT& b) const
    {
        Self pint = integral();

        return pint.value(b) - pint.value(a);
    }

    /* Antiderivative polynomial */
    Self integral(const NT& c = NT(0)) const
    {
        int n = degree() + 1;
        Self result(n);
  
        NT N(n);
        while (n >= 1)
        {
            result.coefs[n] = coefs[n-1] / N;
            n--;
        }
        result.coefs[0] = c;

        //  result.recompute_degree();
        return result;
    }

    /* value at a point */
    NT value(const NT& x) const
    {
        NT result(0);

        for (int i = degree(); i >= 0; --i)
        {
            result *= x;
            result += (*this)[i];
        }

        return result;
    }

    /* compute value and derivative concurrently */
    NT value_and_der(const NT& x, NT& der) const
    {
        NT value = coefs[deg];

        der = 0;
        for (int j = 1; j <= deg; ++j)
        {
            der = der * x + value;
            value = value * x + coefs[deg-j];
        }

        return value;
    }

    Self gcd(const Self& q) const
    {
        if (this->degree() < q.degree()) { return q.gcd(*this); }

        Self h0 = *this, h1 = q, modulo;

        modulo = h0 % h1;
        while ( !modulo.is_zero() )
        {
            h0 = h1;
            h1 = modulo;
            modulo = h0 % h1;
        }

        return h1;
    }
};


template <typename NT>
inline polynomial<NT> gcd(const polynomial<NT>& p,
			  const polynomial<NT>& q)
{
    return p.gcd(q);
}



template <typename NT>
inline polynomial<NT> polymat_det3(polynomial<NT> m[3][3])
{
    polynomial<NT> result;
 
    result  = m[0][2] * ( m[1][0] * m[2][1] - m[1][1] * m[2][0] );
    result -= m[1][2] * ( m[0][0] * m[2][1] - m[0][1] * m[2][0] );
    result += m[2][2] * ( m[0][0] * m[1][1] - m[0][1] * m[1][0] );

    return result;
}



template <typename NT>
inline std::ostream& operator<<(std::ostream &os, const polynomial<NT>& p)
{
    char var = 'x';
    int i, prec;

    prec = os.precision();
    os.precision(polynomial<NT>::IOS_PREC);

    if (!p.degree())
    {
        os << p[p.degree()];
        os.precision(prec);
        return os;
    }

    if (p.degree() > 1)
    {
        if (p[p.degree()] != 1 && p[p.degree()] != -1)
        {
            os << p[p.degree()] << " " << var << "^" << p.degree();
        }
    }
    else
    {
        os << ((p[p.degree()] >= 0) ? "": "-") << var << "^" << p.degree();
    }

    for (i = p.degree() - 1; i > 1; i--)
    {
        if (p[i])
        {
            if (p[i] != 1 && p[i] != -1)
            {
                os << ((p[i] >= 0) ? "+": "") << p[i] << " " << var << "^" << i;
            }
            else
            {
                os << ((p[i] >= 0) ? "+": "-") << var << "^" << i;
            }
        }
    }

    if (p[1])
    {
        if (p[1] != 1 && p[1] != -1)
        {
            if (p.degree() > 1)
            {
                os << ((p[1] >= 0) ? "+": "");
            }
            os << p[1] << " " << var;
        }
        else
        {
            if (p.degree() > 1)
            {
                os << ((p[1] >= 0) ? "+": "-");
            }
            else
            {
                os << ((p[1] >= 0) ? "": "-");
            }
            os << var;
        }
    }
    if (p[0])
    {
        os << ((p[0] >= 0) ? "+": "") << p[0];
    }
    os.precision(prec);
    return os;
}


}} // namespace util::math

} // namespace shape_preserving


#endif // SHAPE_PRESERVING_UTIL_MATH_POLYNOMIAL_HPP
