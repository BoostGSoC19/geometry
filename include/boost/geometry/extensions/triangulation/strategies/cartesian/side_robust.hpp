// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2019 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_SIDE_ROBUST_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_SIDE_ROBUST_HPP

#include <boost/geometry/util/select_most_precise.hpp>
#include <boost/geometry/extensions/triangulation/strategies/cartesian/detail/precise_math.hpp>

namespace boost { namespace geometry
{

namespace strategy { namespace side
{

/*!
\brief Adaptive precision predicate to check at which side of a segment a point lies:
    left of segment (>0), right of segment (< 0), on segment (0).
\ingroup strategies
\tparam CalculationType \tparam_calculation (numeric_limits<ct>::epsilon() and numeric_limits<ct>::digits must be supported for calculation type ct)
\tparam robustness Number that determines maximum precision. Values from 0 to 2 may make the calculation terminate faster for inputs that may require higher precision to ensure correctness.
\details This predicate determines at which side of a segment a point lies using an algorithm that is adapted from orient2d as described in "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates" by Jonathan Richard Shewchuk ( https://dl.acm.org/citation.cfm?doid=237218.237337 ). More information and copies of the paper can also be found at https://www.cs.cmu.edu/~quake/robust.html . It is designed to be adaptive in the sense that it should be fast for inputs that lead to correct results with plain float operations but robust for inputs that require higher precision arithmetics.
 */
template
<
    typename CalculationType = void,
    int robustness = 3
>
struct side_robust
{
public:
    //! \brief Computes double the signed area of the CCW triangle p1, p2, p
    template
    <
        typename CoordinateType,
        typename PromotedType,
        typename P1,
        typename P2,
        typename P
    >
    static inline PromotedType side_value(P1 const& p1, P2 const& p2,
        P const& p)
    {
        std::array<PromotedType, 2> pa {{ get<0>(p1), get<1>(p1) }};
        std::array<PromotedType, 2> pb {{ get<0>(p2), get<1>(p2) }};
        std::array<PromotedType, 2> pc {{ get<0>(p), get<1>(p) }};
        return ::boost::geometry::detail::precise_math::orient2d
            <PromotedType, robustness>(pa, pb, pc);
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template
    <
        typename P1,
        typename P2,
        typename P
    >
    static inline int apply(P1 const& p1, P2 const& p2, P const& p)
    {
        typedef typename coordinate_type<P1>::type coordinate_type1;
        typedef typename coordinate_type<P2>::type coordinate_type2;
        typedef typename coordinate_type<P>::type coordinate_type3;

        typedef typename boost::mpl::if_c
            <
                boost::is_void<CalculationType>::type::value,
                typename select_most_precise
                    <
                        typename select_most_precise
                            <
                                coordinate_type1, coordinate_type2
                            >::type,
                        coordinate_type3
                    >::type,
                CalculationType
            >::type coordinate_type;
        typedef typename select_most_precise
            <
                coordinate_type,
                double
            >::type promoted_type;


        promoted_type sv =
            side_value<coordinate_type, promoted_type>(p1, p2, p);
        return sv > 0 ? 1
            : sv < 0 ? -1
            : 0;
    }
#endif

};

}} // namespace strategy::side

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_SIDE_ROBUST_HPP
