// Shape-preserving interpolation library

// Copyright (c) 2014 Menelaos Karavelas, Heraklion, Greece.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SHAPE_PRESERVING_UTIL_MATH_STANDARD_SEQUENCE_HPP
#define SHAPE_PRESERVING_UTIL_MATH_STANDARD_SEQUENCE_HPP

#include <cassert>
#include <cmath>

#include <shape_preserving/until/math/polynomial.hpp>

namespace shape_preserving
{

namespace util { namespace math
{


template <typename NT_t>
class standard_sequence
{
public:
  typedef NT_t                               NT;
  typedef polynomial<NT>                     polynomial_type;
  typedef standard_sequence<NT>              Self;

private:
    unsigned int nsize;
    polynomial_type p[polynomial::max_degree];

private:
    int sign_variations(const NT& x) const
    {
        NT values[polynomial::max_degree + 1];

        values[0] = p[0].value_and_der(x, values[1]);
        for(unsigned int k = 2; k < nsize; k++) {
            values[k] = p[k].value(x);
        }

        unsigned int count = 0, i = 0, j;
        while (i < nsize && values[i] == 0) { ++i; }
        if (i == nsize) { return 0; }
        while (i < nsize - 1) {
            j = i + 1;
            while (j < nsize && values[j] == 0) { ++j; }
            if (j < nsize && values[i] * values[j] < 0) { ++count; }
            i = j;
        }
        return count;
    }

    void compute_sequence(const polynomial_type& pp)
    {
        assert(pp.degree() < polynomial::max_degree);

        p[0] = pp;
        nsize = 1;

        if (pp.degree() == 0) { return; }
        p[1] = pp.derivative();
        ++nsize;

        polynomial_type modulo = -(p[0] % p[1]);

        while ( !modulo.is_zero() ) {
            p[nsize] = modulo;
            modulo = -(p[nsize - 1] % p[nsize]);
            ++nsize;
        }
    }

public:
    standard_sequence() { nsize = 0; }
    standard_sequence(const polynomial_type &pp)
    {
        compute_sequence(pp);
    }

    standard_sequence(const Self& seq)
    {
        nsize = seq.nsize;
        for (unsigned int i = 0; i < nsize; ++i) {
            p[i] = seq.p[i];
        }
    }

    inline void set_polynomial(const polynomial_type& pp)
    {
        compute_sequence(pp);
    }

    inline unsigned int size() { return nsize; }

    inline polynomial_type operator[](unsigned int i) const { return p[i]; }

    int number_of_real_roots(const NT& a, const NT& b,
                             unsigned int max_nr) const
    {
        if (max_nr == 0) { return 0; } // THIS SHOULD NOT BE HERE
        NT pa = p[0].value(a), pb = p[0].value(b);
        assert(pa != 0 && pb != 0);
        if (max_nr == 1 && p[nsize - 1].degree() == 0) {
            if (pa * pb < 0) { return 1; }
            else { return 0; }
        } else if (max_nr > 1) {
            int sv1 = sign_variations(a);
            if (!sv1) { return 0; }
            int sv = sv1 - sign_variations(b);
            return ((sv >= 0) ? sv : 0);
        }
        return 0;
    }
};


}} // namespace util::math

} // namespace shape_preserving



#endif // SHAPE_PRESERVING_UTIL_MATH_STANDARD_SEQUENCE_HPP
