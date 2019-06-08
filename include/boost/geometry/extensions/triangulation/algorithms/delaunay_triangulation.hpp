#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_ALGORITHMS_DELAUNAY_TRIANGULATION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_ALGORITHMS_DELAUNAY_TRIANGULATION_HPP

#include <vector>
#include <utility>
#include <tuple>
#include <iterator>
#include <algorithm>
#include <boost/geometry/extensions/triangulation/geometries/triangulation.hpp>
#include <boost/geometry/extensions/triangulation/strategies/cartesian/in_circle_by_determinant.hpp>
#include <boost/geometry/strategies/cartesian/side_by_triangle.hpp>

namespace boost { namespace geometry
{

#ifndef DOXYGEN_NO_DETAIL

namespace detail { namespace delaunay_triangulation
{

template
<
    typename Area,
    typename Point
>
Area triangle_area(const Point& p1, const Point& p2, const Point& p3) 
{
    return determinant<Area>(
        get<0>(p2)-get<0>(p1), get<1>(p2)-get<1>(p1),
        get<0>(p3)-get<0>(p1), get<1>(p3)-get<1>(p1))/2;
}

template<typename Point>
typename default_comparable_distance_result<Point, Point>::type comparable_circumcircle_diameter(const Point& p1, const Point& p2, const Point& p3)
{
    typedef typename default_comparable_distance_result<Point, Point>::type dfcr;
    dfcr comp_area = triangle_area<dfcr>(p1,p2,p3);
    return comparable_distance(p1,p2) * comparable_distance(p1,p3) * comparable_distance(p2,p3)
        / (comp_area * comp_area);
}

template
<
    typename Point
>
Point circumcircle_center(const Point& p1, const Point& p2, const Point& p3)
{
    auto ax = get<0>(p1);
    auto ay = get<1>(p1);
    auto bx = get<0>(p2);
    auto by = get<1>(p2);
    auto cx = get<0>(p3);
    auto cy = get<1>(p3);
    auto d = 2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
    auto x = ((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d;
    auto y = ((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d;
    return boost::geometry::make<Point>(x,y);
}

template<typename PointContainer, typename Triangulation, typename SideStrategy, typename InCircleStrategy, typename CalculationType=double> // TODO: void CT
void delaunay_triangulation(PointContainer const & in, Triangulation& out, bool legalize)
{
    typedef typename PointContainer::value_type point_type;
    typedef typename default_distance_result<point_type, point_type>::type distance_type;
    typedef CalculationType ct;
    std::vector<std::tuple<point_type, ct, ct>> points;
    points.reserve(std::size(in));
    //Step 1
    points.emplace_back(*std::begin(in),ct(0), ct(0));
    //Step 2 & 3
    std::transform(std::begin(in) + 1, std::end(in), std::back_inserter(points),
        [&points](point_type const& p) {
            return std::tuple<point_type, ct, ct>(p, 
            boost::geometry::comparable_distance(p, std::get<0>(points[0])),
            0);
        });
    std::sort(std::begin(points)+1, std::end(points), 
            [](std::tuple<point_type, ct, ct> const& p0, std::tuple<point_type, ct, ct> const& p1) 
            { return std::get<1>(p0) < std::get<1>(p1); });

    //Step 4
    {
        distance_type min_circumdiameter = std::numeric_limits<distance_type>::max();
        std::size_t min_index = 2;
        for( std::size_t i = 2 ; std::get<1>(points[i]) < min_circumdiameter && i < points.size() ; ++i)
        {
            distance_type diam = comparable_circumcircle_diameter(
                std::get<0>(points[0]), std::get<0>(points[1]), std::get<0>(points[i]));
            if(diam < min_circumdiameter) {
                min_index = i;
                min_circumdiameter = diam;
            }
        }
        std::swap(points[2], points[min_index]);
    }
    //Step 5
    if(SideStrategy::apply(std::get<0>(points[0]), std::get<0>(points[1]), std::get<0>(points[2])) < 0) {
        std::swap(points[1], points[2]);
    }
    const auto PI = 3.14159265358979323846;
    //Step 6
    
    const auto p1 = out.add_vertex(std::get<0>(points[0]));                             
    const auto p2 = out.add_vertex(std::get<0>(points[1]));
    const auto p3 = out.add_vertex(std::get<0>(points[2]));
    auto seed_face = out.add_isolated_face(p1, p2, p3);
    point_type cen;
    set<0>(cen, get<0>(std::get<0>(points[0]))/3 + get<0>(std::get<0>(points[1]))/3 + get<0>(std::get<0>(points[2]))/3);
    set<1>(cen, get<1>(std::get<0>(points[0]))/3 + get<1>(std::get<0>(points[1]))/3 + get<1>(std::get<0>(points[2]))/3);

    point_type C = circumcircle_center( std::get<0>(points[0]), std::get<0>(points[1]), std::get<0>(points[2]));
    {
        std::for_each(std::begin(points), std::end(points),
            [&C, &cen, &PI](std::tuple<point_type, ct, ct>& p){ 
                std::get<1>(p) = comparable_distance(std::get<0>(p), C);
                std::get<2>(p) = std::atan2( get<1>(std::get<0>(p)) - get<1>(cen), get<0>(std::get<0>(p)) - get<0>(cen) ) + PI;
            });
        std::sort(std::begin(points)+3, std::end(points),
            [](std::tuple<point_type, ct, ct> const& p0, std::tuple<point_type, ct, ct> const& p1)
            { return std::get<1>(p0) < std::get<1>(p1); });
    }
    //Step 7 & 8
    {
        auto e1 = out.face_edge(seed_face);
        auto e2 = out.next(e1);
        auto e3 = out.next(e2);
        std::vector<std::pair<typename Triangulation::halfedge_index, ct>> convex_hull { 
            std::make_pair(e1, ct(std::get<2>(points[1]))), 
            std::make_pair(e2, ct(std::get<2>(points[2]))), 
            std::make_pair(e3, ct(std::get<2>(points[0]))) };
        for(int i=3; i<points.size(); ++i)
        {
            auto new_vertex = out.add_vertex(std::get<0>(points[i]));
            auto const& p = out.vertex(new_vertex);
            auto is_visible = [&out, &p](std::pair<typename Triangulation::halfedge_index, ct> const& bep)
                {
                    const auto& be = std::get<0>(bep);
                    const auto s = out.face_segment(be);
                    const auto p1 = s.first;
                    const auto p2 = s.second;
                    bool result = SideStrategy::apply(
                        p1,
                        p2,     
                        p)<0;
                    return result;
                };
            const auto& ref_angle = std::get<1>(*convex_hull.begin());
            auto angle = std::get<2>(points[i])-ref_angle;       
            if(angle < 0) angle += 2*PI;
            ct opposite = angle - PI;
            if(opposite < 0) opposite += 2*PI;

            for(int i=0; i < convex_hull.size() ; ++i) {
                auto ba = std::get<1>(convex_hull[i]);
                auto p = out.face_segment(std::get<0>(convex_hull[i])).first;
            }  
            for(int i=0; i < convex_hull.size() ; ++i) {
                auto ba = std::get<1>(convex_hull[i]) - ref_angle;
                if(ba < 0) ba += 2*PI;
            }
            auto pred = [&ref_angle, &PI]
                (std::pair<typename Triangulation::halfedge_index, ct> const& be, ct const& a)
                {
                    auto ba = std::get<1>(be) - ref_angle;
                    if(ba < 0) ba += 2*PI;
                    return ba < a;
                };
            auto vis_edge = std::lower_bound(convex_hull.begin(), convex_hull.end(), angle, pred);
            if( vis_edge == convex_hull.begin() ) vis_edge = convex_hull.end()-1; else --vis_edge;
            auto invis_edge = std::lower_bound(convex_hull.begin(), convex_hull.end(), opposite, pred);
            if( invis_edge == convex_hull.begin() ) invis_edge = convex_hull.end()-1; else --invis_edge;

            auto vis_small = [&is_visible] 
                (std::pair<typename Triangulation::halfedge_index, ct> const& be, int const&)
                {
                    if(is_visible(be)) return true;
                    else return false;
                };
            auto vis_large = [&is_visible]
                (std::pair<typename Triangulation::halfedge_index, ct> const& be, int const&)
                {
                    if(is_visible(be)) return false;
                    else return true; 
                };
            typedef typename std::vector<std::pair<typename Triangulation::halfedge_index, ct>>::iterator ch_iterator;
            ch_iterator first_visible, last_visible, fv2 = convex_hull.end();
            if(vis_edge < invis_edge) {
                if(is_visible(*convex_hull.begin())) {
                    first_visible = convex_hull.begin();
                    last_visible = std::lower_bound(vis_edge, invis_edge, 0, vis_small);
                    fv2 = std::lower_bound(invis_edge, convex_hull.end(), 0, vis_large);
                } else {
                    first_visible = std::lower_bound(convex_hull.begin(), vis_edge, 0, vis_large);
                    last_visible = std::lower_bound(vis_edge, convex_hull.end(), 0, vis_small);
                }
            } else {
                if(is_visible(*convex_hull.begin())) {
                    first_visible = convex_hull.begin();
                    last_visible = std::lower_bound(convex_hull.begin()+1, invis_edge, 0, vis_small);
                    fv2 = std::lower_bound(invis_edge, vis_edge, 0, vis_large);
                } else {
                    first_visible = std::lower_bound(invis_edge, vis_edge, 0, vis_large);
                    last_visible = std::lower_bound(vis_edge, convex_hull.end(), 0, vis_small);
                }
            }
            auto first_visible_angle = std::get<1>(*first_visible);
            const bool begin_visible = first_visible == std::begin(convex_hull);
            typename Triangulation::face_index prev = -1, fnf2 = -1;
            bool looped = false;
            ct fv2_angle;
            if(begin_visible && is_visible(convex_hull.back()))
            {
                looped = true;
                fv2_angle = std::get<1>(*fv2);
                fnf2 = out.add_face_on_boundary(std::get<0>(*fv2), new_vertex);
                prev = fnf2;
                for(auto it = fv2 + 1; it != std::end(convex_hull); ++it)
                {
                    auto next = out.add_face_on_boundary(std::get<0>(*it), i);
                    out.connect(out.next(out.face_edge(next)), out.prev(out.face_edge(prev)));
                    prev = next;
                }
                convex_hull.erase(fv2, std::end(convex_hull));
            }
            const auto fnf = out.add_face_on_boundary(std::get<0>(*first_visible), new_vertex);
            if(looped) {
                out.connect(out.next(out.face_edge(fnf)), out.prev(out.face_edge(prev)));
            }
            prev = fnf;
            for(auto it = first_visible + 1; it != last_visible; ++it)
            {
                auto next = out.add_face_on_boundary(std::get<0>(*it), i);
                out.connect(out.next(out.face_edge(next)), out.prev(out.face_edge(prev)));
                prev = next;
            }
            auto ip = convex_hull.erase(first_visible, last_visible);
            ip = convex_hull.insert(ip, std::make_pair(out.prev(out.face_edge(prev)), ct(std::get<2>(points[i]))));
            if(looped)
                convex_hull.insert(ip, std::make_pair(out.next(out.face_edge(fnf2)), ct(fv2_angle)));
            else
                convex_hull.insert(ip, std::make_pair(out.next(out.face_edge(fnf)), ct(first_visible_angle)));
        }
    }
    //Step 9
    if(legalize)
    {
        std::vector<typename Triangulation::halfedge_index> L;
        for(std::size_t i = 0; i < out.faces(); ++i)
        {
            for(unsigned short j = 0; j<2 ; ++j)
                if(i > out.neighbour(i, j)) {
                    L.push_back(typename Triangulation::halfedge_index(i,j));
                }
        }
        auto edge_legal = [&out](typename Triangulation::halfedge_index const& e)
            {
                auto const& p1 = out.face_vertex(e.m_f, 0);
                auto const& p2 = out.face_vertex(e.m_f, 1);
                auto const& p3 = out.face_vertex(e.m_f, 2);
                auto const& p = out.face_vertex(
                        out.neighbour(e.m_f, e.m_v),
                        out.opposite(e.m_f, e.m_v)
                    );
                return InCircleStrategy::apply(p1, p2, p3, p) <= 0;
            };
        while( !std::empty(L) )
        {
            auto e = L.back();
            L.pop_back();
            auto const& f1 = e.m_f;
            auto const& v1 = e.m_v;
            auto const f2 = out.neighbour(f1, v1);
            if(f2 == Triangulation::invalid)
                continue;
            auto const v2 = out.opposite(f1, v1);
            if( !edge_legal(e) ) {
                L.emplace_back(f1, v1);
                L.emplace_back(f2, v2);
                L.emplace_back(f1, v1 == 0 ? 2 : v1 - 1);
                L.emplace_back(f2, v2 == 0 ? 2 : v2 - 1);
                out.flip(e);
            }
        }
    }
}

}} // namespace detail::delaunay_triangulation

template<typename PointContainer, typename Triangulation, typename SideStrategy = strategy::side::side_by_triangle<>, typename InCircleStrategy = strategy::in_circle::fast_in_circle<>>
void delaunay_triangulation(PointContainer const & in, Triangulation& out, bool legalize = true)
{
    detail::delaunay_triangulation::delaunay_triangulation<PointContainer, Triangulation, SideStrategy, InCircleStrategy>(in, out, legalize);
}

#endif // DOXYGEN_NO_DETAIL

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_ALGORITHMS_DELAUNAY_TRIANGULATION_HPP
