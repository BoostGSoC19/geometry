#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_SIDE_ROBUST_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_SIDE_ROBUST_HPP
#include<boost/geometry/extensions/triangulation/strategies/cartesian/detail/precise_math.hpp>

namespace boost { namespace geometry
{ 

namespace strategy { namespace side
{

template <typename CalculationType = double>
struct side_robust
{
public:
    template <typename P1, typename P2, typename P3>
    static inline CalculationType side_value(P1 const& p1, P2 const& p2, P3 const& p3)
    {
        std::array<CalculationType, 2> pa {{ boost::geometry::get<0>(p1), boost::geometry::get<1>(p1) }};
        std::array<CalculationType, 2> pb {{ boost::geometry::get<0>(p2), boost::geometry::get<1>(p2) }};
        std::array<CalculationType, 2> pc {{ boost::geometry::get<0>(p3), boost::geometry::get<1>(p3) }};
        return boost::geometry::detail::precise_math::orient2d<CalculationType>(pa, pb, pc);
    }

    template <typename P1, typename P2, typename P3>
    static inline int apply(P1 const& p1, P2 const& p2, P3 const& p3)
    {
        CalculationType sv = side_value(p1, p2, p3);
        return sv > 0 ? 1
            : sv < 0 ? -1
            : 0;
    }

};

}} // namespace strategy::side

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_SIDE_ROBUST_HPP

