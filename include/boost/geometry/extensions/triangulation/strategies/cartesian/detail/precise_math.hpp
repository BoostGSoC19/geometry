#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_DETAIL_PRECISE_MATH_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_DETAIL_PRECISE_MATH_HPP
#include<numeric>
#include<array>
// The following code is based on 
// "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates" by 
// Richard Shewchuk, J. Discrete Comput Geom (1997) 18: 305. https://doi.org/10.1007/PL00009321

namespace boost { namespace geometry 
{ 

namespace detail { namespace precise_math
{

// See Theorem 6, page 6
template
<
    typename RealNumber
>
inline std::array<RealNumber, 2> fast_two_sum(RealNumber const a, RealNumber const b)
{
    RealNumber x = a + b;
    RealNumber b_virtual = x - a;
    return {{x, b - b_virtual}};
}

// See Theorem 7, page 7 - 8
template
<
    typename RealNumber
>
inline std::array<RealNumber, 2> two_sum(RealNumber const a, RealNumber const b)
{
    RealNumber x = a + b;
    RealNumber b_virtual = x - a;
    RealNumber a_virtual = x - b_virtual;
    RealNumber b_roundoff = b - b_virtual;
    RealNumber a_roundoff = a - a_virtual;
    RealNumber y = a_roundoff + b_roundoff;
    return {{ x,  y }};
}

// See bottom of page 8
template
<
    typename RealNumber
>
inline RealNumber two_diff_tail(RealNumber const a, RealNumber const b, RealNumber const x)
{
    RealNumber b_virtual = a - x;
    RealNumber a_virtual = x + b_virtual;
    RealNumber b_roundoff = b_virtual - b;
    RealNumber a_roundoff = a - a_virtual;
    return a_roundoff + b_roundoff;
}

// see bottom of page 8
template
<
    typename RealNumber
>
inline std::array<RealNumber, 2> two_diff(RealNumber const a, RealNumber const b)
{
    RealNumber x = a - b;
    RealNumber y = two_diff_tail(a, b, x);
    return {{ x, y }};
}

// constexpr power-method, helper for splitter
template
<
    typename RealNumber
>
constexpr RealNumber int_pow(RealNumber const base, int exp, RealNumber out = 1.0)
{
    return exp < 1 ? out : int_pow<RealNumber>(base*base, exp/2, (exp % 2) ? out*base : out);
}

// consexpr method to compute 2^s + 1 as in Theorem 17, page 18
template
<
    typename RealNumber
>
constexpr RealNumber splitter()
{
    return int_pow<RealNumber>(2.0, (std::numeric_limits<RealNumber>::digits+1)/2) + 1;
}

// see theorem 17, page 18
template
<
    typename RealNumber
>
inline std::array<RealNumber, 2> split(RealNumber const a) {
    RealNumber c = splitter<RealNumber>() * a;
    RealNumber a_big = c - a;
    RealNumber a_hi = c - a_big;
    return {{ a_hi, a - a_hi }};
}

// see theorem 18, page 19
template
<
    typename RealNumber
>
inline RealNumber two_product_tail(RealNumber const a, RealNumber const b, RealNumber const x)
{
    std::array<RealNumber, 2> a_expansion = split(a);
    std::array<RealNumber, 2> b_expansion = split(b);
    RealNumber err1 = x - (a_expansion[0] * b_expansion[0]);
    RealNumber err2 = err1 - (a_expansion[1] * b_expansion[0]);
    RealNumber err3 = err2 - (a_expansion[0] * b_expansion[1]);
    return (a_expansion[1] * b_expansion[1]) - err3;
}

// see theorem 18, page 19
template
<
    typename RealNumber
>
inline std::array<RealNumber, 2> two_product(RealNumber const a, RealNumber const b)
{
    RealNumber x = a * b;
    RealNumber y = two_product_tail(a, b, x);
    return {{ x , y }};
}

// see theorem 12, figure 7, page 11 - 12, 
// this is the 2 by 2 case for the corresponding diff-method
// note that this method takes input in descending order of magnitude and
// returns components in ascending order of magnitude
template
<
    typename RealNumber
>
inline std::array<RealNumber, 4> two_two_expansion_diff(std::array<RealNumber, 2> const a,
        std::array<RealNumber, 2> const b)
{
    std::array<RealNumber, 4> h;
    std::array<RealNumber, 2> Qh = two_diff(a[1], b[1]);
    h[0] = Qh[1];
    Qh = two_sum( a[0], Qh[0] );
    RealNumber _j = Qh[0];
    Qh = two_diff(Qh[1], b[0]);
    h[1] = Qh[1];
    Qh = two_sum( _j, Qh[0] );
    h[2] = Qh[1];
    h[3] = Qh[0];
    return h;
}

// see theorem 13, figure 8. This implementation uses zero elimination as
// suggested on page 17, second to last paragraph. Returns the number of
// non-zero components in the result and writes the result to h.
// the merger into a single sequence g is done implicitly
template
<
    typename RealNumber,
    std::size_t InSize1,
    std::size_t InSize2,
    std::size_t OutSize
>
inline int fast_expansion_sum_zeroelim(std::array<RealNumber, InSize1> const& e,
        std::array<RealNumber, InSize2> const& f, std::array<RealNumber, OutSize> & h,
        int m = InSize1, int n = InSize2)
{
    std::array<RealNumber, 2> Qh;
    int i_e = 0, i_f = 0, i_h = 0;
    if (std::abs(f[0]) > std::abs(e[0])) {
        Qh[0] = e[i_e++];
    } else {
        Qh[0] = f[i_f++];
    }
    i_h = 0;
    if ((i_e < m) && (i_f < n)) {
        if (std::abs(f[i_f]) > std::abs(e[i_e])) {
            Qh = fast_two_sum(e[i_e++], Qh[0]);
        } else {
            Qh = fast_two_sum(f[i_f++], Qh[0]);
        }
        if (Qh[1] != 0.0) {
            h[i_h++] = Qh[1];
        }
        while ((i_e < m) && (i_f < n)) {
            if (std::abs(f[i_f]) > std::abs(e[i_e])) {
                Qh = two_sum(Qh[0], e[i_e++]);
            } else {
                Qh = two_sum(Qh[0], f[i_f++]);
            }
            if (Qh[1] != 0.0) {
                h[i_h++] = Qh[1];
            }
        }
    }
    while (i_e < m) {
        Qh = two_sum(Qh[0], e[i_e++]);
        if (Qh[1] != 0.0) {
            h[i_h++] = Qh[1];
        }
    }
    while (i_f < n) {
        Qh = two_sum(Qh[0], f[i_f++]);
        if (Qh[1] != 0.0) {
            h[i_h++] = Qh[1];
        }
    }
    if ((Qh[0] != 0.0) || (i_h == 0)) {
        h[i_h++] = Qh[0];
    }
    return i_h;
}

// see theorem 19, figure 13, page 20 - 21. This implementation uses zero
// elimination as suggested on page 17, second to last paragraph. Returns the
// number of non-zero components in the result and writes the result to h.
template
<
    typename RealNumber,
    std::size_t InSize
>
inline int scale_expansion_zeroelim(std::array<RealNumber, InSize> const& e,
        RealNumber const b, std::array<RealNumber, 2 * InSize> & h, int e_non_zeros = InSize)
{
    std::array<RealNumber, 2> Qh = two_product(e[0], b);
    int i_h = 0;
    if (Qh[1] != 0) {
        h[i_h++] = Qh[1];
    }
    for (int i_e = 1; i_e < e_non_zeros; i_e++) {
        std::array<RealNumber, 2> Tt = two_product(e[i_e], b);
        Qh = two_sum(Qh[0], Tt[1]);
        if (Qh[1] != 0) {
            h[i_h++] = Qh[1];
        }
        Qh = fast_two_sum(Tt[0], Qh[0]);
        if (Qh[1] != 0) {
            h[i_h++] = Qh[1];
        }
    }
    if ((Qh[0] != 0.0) || (i_h == 0)) {
        h[i_h++] = Qh[0];
    }
    return i_h;
}

// see page 38, Figure 21 for the calculations, notation follows the notation in the figure.
template
<
    typename RealNumber
>
inline RealNumber orient2d(std::array<RealNumber, 2> const& p1,
        std::array<RealNumber, 2> const& p2, std::array<RealNumber, 2> const& p3)
{
    std::array<RealNumber, 2> t1, t2, t3, t4;
    t1[0] = p1[0] - p3[0];
    t2[0] = p2[1] - p3[1];
    t3[0] = p1[1] - p3[1];
    t4[0] = p2[0] - p3[0];
    std::array<RealNumber, 2> t5_01, t6_01;
    t5_01[0] = t1[0] * t2[0];
    t6_01[0] = t3[0] * t4[0];
    RealNumber det = t5_01[0] - t6_01[0];

    if ( (t5_01[0] > 0 && t6_01[0] <= 0) || (t5_01[0] < 0 && t6_01[0] >= 0) ) {
        //if diagonal and antidiagonal have different sign, the sign of det is obvious
        return det;
    }
    RealNumber const magnitude = std::abs(t5_01[0]) + std::abs(t6_01[0]);

    // see p.39, mind the different definition of epsilon for error bound
    RealNumber const A_relative_bound = (1.5 + 4 * std::numeric_limits<RealNumber>::epsilon())
        * std::numeric_limits<RealNumber>::epsilon();
    RealNumber absolute_bound = A_relative_bound * magnitude;
    if ( std::abs(det) >= absolute_bound ) {
        return det; //A estimate
    }

    t5_01[1] = two_product_tail(t1[0], t2[0], t5_01[0]);
    t6_01[1] = two_product_tail(t3[0], t4[0], t6_01[0]);
    std::array<RealNumber, 4> tA_03 = two_two_expansion_diff(t5_01, t6_01);
    det = std::accumulate(tA_03.begin(), tA_03.end(), static_cast<RealNumber>(0));
    // see p.39, mind the different definition of epsilon for error bound
    RealNumber B_relative_bound = (1 + 3 * std::numeric_limits<RealNumber>::epsilon())
        * std::numeric_limits<RealNumber>::epsilon();
    absolute_bound = B_relative_bound * magnitude;
    if (std::abs(det) >= absolute_bound) {
        return det; //B estimate
    }
    t1[1] = two_diff_tail(p1[0], p3[0], t1[0]);
    t2[1] = two_diff_tail(p2[1], p3[1], t2[0]);
    t3[1] = two_diff_tail(p1[1], p3[1], t3[0]);
    t4[1] = two_diff_tail(p2[0], p3[0], t4[0]);

    if ((t1[1] == 0.0) && (t3[1] == 0.0) && (t2[1] == 0.0) && (t4[1] == 0.0)) {
        return det; //If all tails are zero, there is noething else to compute
    }
    RealNumber sub_bound = (1.5 + 2 * std::numeric_limits<RealNumber>::epsilon())
        * std::numeric_limits<RealNumber>::epsilon();
    // see p.39, mind the different definition of epsilon for error bound
    RealNumber C_relative_bound = (2.25 + 8 * std::numeric_limits<RealNumber>::epsilon()) 
        * std::numeric_limits<RealNumber>::epsilon() * std::numeric_limits<RealNumber>::epsilon();
    absolute_bound = C_relative_bound * magnitude + sub_bound * std::abs(det);
    det += (t1[0] * t2[1] + t2[0] * t1[1]) - (t3[0] * t4[1] + t4[0] * t3[1]);
    if (std::abs(det) >= absolute_bound) {
        return det; //C estimate
    }
    std::array<RealNumber, 8> D_left;
    int D_left_nz;
    {
        std::array<RealNumber, 2> t5_23 = two_product(t1[1], t2[0]);
        std::array<RealNumber, 2> t6_23 = two_product(t3[1], t4[0]);
        std::array<RealNumber, 4> tA_47 = two_two_expansion_diff(t5_23, t6_23);
        D_left_nz = fast_expansion_sum_zeroelim(tA_03, tA_47, D_left);
    }
    std::array<RealNumber, 8> D_right;
    int D_right_nz;
    {
        std::array<RealNumber, 2> t5_45 = two_product(t1[0], t2[1]);
        std::array<RealNumber, 2> t6_45 = two_product(t3[0], t4[1]);
        std::array<RealNumber, 4> tA_8_11 = two_two_expansion_diff(t5_45, t6_45);
        std::array<RealNumber, 2> t5_67 = two_product(t1[1], t2[1]);
        std::array<RealNumber, 2> t6_67 = two_product(t3[1], t4[1]);
        std::array<RealNumber, 4> tA_12_15 = two_two_expansion_diff(t5_67, t6_67);
        D_right_nz = fast_expansion_sum_zeroelim(tA_8_11, tA_12_15, D_right);
    }
    std::array<RealNumber, 16> D;
    int D_nz = fast_expansion_sum_zeroelim(D_left, D_right, D, D_left_nz, D_right_nz);
    // only return component of highest magnitude because we mostly care about the sign.
    return(D[D_nz - 1]);
}

}} // namespace detail::precise_math

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_DETAIL_PRECISE_MATH_HPP
