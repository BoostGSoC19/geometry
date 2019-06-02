// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2019 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_RANDOM_UNIFORM_POINT_DISTRIBUTION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_RANDOM_UNIFORM_POINT_DISTRIBUTION_HPP

#include <boost/geometry/extensions/random/dispatch/uniform_point_distribution.hpp>

namespace boost { namespace geometry { namespace random
{

    template<typename DomainGeometry, typename Point = typename geometry::point_type<DomainGeometry>::type>
    dispatch::uniform_point_distribution<DomainGeometry, Point>
    uniform_point_distribution(DomainGeometry const& domain)
    {
        return dispatch::uniform_point_distribution<DomainGeometry, Point>(domain);
    }

}}} // namespace boost::geometry::random

#endif // BOOST_GEOMETRY_EXTENSIONS_RANDOM_UNIFORM_POINT_DISTRIBUTION_HPP
