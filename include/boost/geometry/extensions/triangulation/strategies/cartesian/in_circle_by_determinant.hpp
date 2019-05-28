#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_POINT_IN_CIRCLE_BY_DETERMINANT_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_POINT_IN_CIRCLE_BY_DETERMINANT_HPP

namespace boost { namespace geometry 
{ 

namespace detail
{
    
template <typename ct>
ct determinant(
        const ct& v11, const ct& v12, const ct& v13, 
        const ct& v21, const ct& v22, const ct& v23, 
        const ct& v31, const ct& v32, const ct& v33)
{
    return v11*v22*v33 + v12*v23*v31 + v13*v21*v32 - v13*v22*v31 - v12*v21*v33 - v11*v23*v32;
}

} //namespace detail

namespace strategy
{

namespace in_circle
{

template <typename CalculationType = void>
class fast_in_circle
{
public:
    template <typename P1, typename P2, typename P3, typename P>
    static inline int apply(P1 const& p1, P2 const& p2, P3 const& p3, P const& p)
    {
        typedef typename coordinate_type<P1>::type coordinate_type1;
        typedef typename coordinate_type<P2>::type coordinate_type2;
        typedef typename coordinate_type<P3>::type coordinate_type3;
        typedef typename coordinate_type<P>::type coordinate_type4;

        typedef typename boost::mpl::if_c 
            <
                boost::is_void<CalculationType>::type::value,
                typename select_most_precise
                    <
                        typename select_most_precise
                            <
                                typename select_most_precise
                                    <
                                        coordinate_type1, coordinate_type2
                                    >::type, 
                                coordinate_type3
                            >::type,
                        coordinate_type4
                    >::type,
                CalculationType
            >::type coordinate_type;
        typedef typename select_most_precise
            <
                coordinate_type,
                double
            >::type promoted_type;
        promoted_type const d1x = get<0>(p1) - get<0>(p);
        promoted_type const d1y = get<1>(p1) - get<1>(p);
        promoted_type const d2x = get<0>(p2) - get<0>(p);
        promoted_type const d2y = get<1>(p2) - get<1>(p);
        promoted_type const d3x = get<0>(p3) - get<0>(p);
        promoted_type const d3y = get<1>(p3) - get<1>(p);
        promoted_type in = geometry::detail::determinant<promoted_type>
            (
                d1x, d1y, d1x*d1x + d1y*d1y,
                d2x, d2y, d2x*d2x + d2y*d2y,
                d3x, d3y, d3x*d3x + d3y*d3y
            );
        promoted_type const zero = promoted_type();
        return in > zero ? 1
            : in < zero ? -1
            : 0;
    }

};

} // namespace in_circle

} // namespace strategy

}} // namespace boost::geometry::strategy

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_STRATEGIES_CARTESIAN_POINT_IN_CIRCLE_BY_DETERMINANT_HPP

