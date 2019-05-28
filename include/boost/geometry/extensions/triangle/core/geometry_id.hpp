#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGLE_CORE_GEOMETRY_ID_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGLE_CORE_GEOMETRY_ID_HPP


#include <boost/geometry/core/geometry_id.hpp>
#include <boost/geometry/extensions/triangle/core/tags.hpp>


namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DISPATCH
namespace core_dispatch
{

template <>
struct geometry_id<triangle_tag> : boost::mpl::int_<95> {};

} // namespace core_dispatch
#endif

}} //namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGLE_CORE_GEOMETRY_ID_HPP
