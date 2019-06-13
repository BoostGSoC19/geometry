#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_GEOMETRIES_TRIANGULATION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_GEOMETRIES_TRIANGULATION_HPP

#include <vector>

#include <boost/range.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/indirect_iterator.hpp>

#include <boost/geometry/strategies/cartesian/side_by_triangle.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/views/detail/points_view.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>

#include <boost/geometry/extensions/triangle/triangle.hpp>

namespace boost { namespace geometry
{
template<typename Triangulation>
struct face_range_type {};

template<typename Triangulation>
typename face_range_type<Triangulation>::type face_range(Triangulation const& t);

template<typename Triangulation>
struct triangulation_face {};

namespace model
{

template
<
    typename Point,
    template<typename, typename> class VertexContainer = std::vector,
    template<typename, typename> class FaceContainer = std::vector,
    template<typename> class VertexAllocator = std::allocator,
    template<typename> class FaceAllocator = std::allocator
>
class triangulation;

template<typename Triangulation>
struct vertex_ref
{
    typename Triangulation::point_type m_p;
    typename Triangulation::face_iterator m_f;
};

template<typename Triangulation>
struct face_ref
{
    typedef typename Triangulation::point_type point_type;
    typedef typename Triangulation::vertex_iterator vertex_iterator;
    typedef typename Triangulation::face_iterator face_iterator;
    std::array<vertex_iterator, 3> m_v;
    std::array<face_iterator, 3> m_f;
    std::array<unsigned short, 3> m_o;
public:
    typedef boost::indirect_iterator<typename std::array<vertex_iterator, 3>::const_iterator>
        const_iterator;
    typedef boost::indirect_iterator<typename std::array<vertex_iterator, 3>::iterator>
        iterator;

    const_iterator begin() const { return const_iterator(m_v.begin()); }
    const_iterator end() const { return const_iterator(m_v.end()); }
};

template
<
    typename Value,
    template<typename, typename> class Container,
    template<typename> class Allocator
>
struct reserve_if_vector{
    static void apply(Container<Value, Allocator<Value>>& c, std::size_t n) {}
};

template<typename Value, template<typename> class Allocator>
struct reserve_if_vector<Value, std::vector, Allocator>{
    static void apply(std::vector<Value, Allocator<Value>>& c, std::size_t n) { c.reserve(n);}
};

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
class triangulation
{
private:
    BOOST_CONCEPT_ASSERT( (concepts::Point<Point>) );
public:
    typedef Point point_type;
    typedef face_ref<triangulation> face_type;
    typedef vertex_ref<triangulation> vertex_type;
    typedef unsigned short face_vertex_index;
    typedef VertexContainer<vertex_type, VertexAllocator<vertex_type>> vertex_container;
    typedef FaceContainer<face_type, FaceAllocator<face_type>> face_container;
    typedef typename face_container::iterator face_iterator;
    typedef typename face_container::const_iterator const_face_iterator;
    typedef typename vertex_container::iterator vertex_iterator;
    typedef typename vertex_container::const_iterator const_vertex_iterator;
    typedef typename coordinate_type<Point>::type coordinate_type;
    typedef typename model::segment<Point> segment_type;
    struct halfedge_index
    {
        halfedge_index(face_iterator f, face_vertex_index v):m_f(f), m_v(v) {}
        face_iterator m_f;
        face_vertex_index m_v;
    };

    struct fulledge_index
    {
        fulledge_index(
            face_iterator f1, face_vertex_index v1, face_iterator f2, face_vertex_index v2)
            :m_f1(f1), m_f2(f2), m_v1(v1), m_v2(v2) {}
        fulledge_index(halfedge_index e)
            :m_f1(e.m_f), m_f2(e.m_f->m_f[e.m_v]), m_v1(e.m_v), m_v2(e.m_f->m_o[e.m_v]) {}
        face_iterator m_f1, m_f2;
        face_vertex_index m_v1, m_v2;
    };

/*    void debug_print()
    {
        std::cout << "Vertices: \n";
        for(std::size_t i = 0; i < m_vertices.size(); ++i)
            std::cout << "Vertex " << i << ": ( " << get<0>(m_vertices[i].m_p) << " , " << get<1>(m_vertices[i].m_p) 
            << "), touches face: " << m_vertices[i].m_f << "\n";
        std::cout << "Faces: \n";
        for(std::size_t i = 0; i < m_faces.size(); ++i) {
            debug_print_face(i);
        }
        std::cout << "boundary vertex: " 
            << std::distance<const_face_iterator>(m_vertices.begin(), m_boundary_vertex) << "\n";
    }

    void debug_print_face(std::size_t i) {
            const face_ref& f = m_faces[i];
            std::cout << "Face " << i << ":\n";
            for(unsigned short v = 0; v < 3; ++v)
                std::cout << "Vertex " << v << ": " 
                    << std::distance<const_vertex_iterator>(m_vertices.begin(), f.m_v[v])
                    << ", Neighbour: " << std::distance<const_face_iterator>(m_faces.begin(), f.m_f[v]) 
                    << ", Opposite: " << f.m_o[v] << "\n";
    }*/

    const face_iterator invalid = face_iterator();
    triangulation(std::size_t points = 3)
    {
        reserve_if_vector<vertex_type, VertexContainer, VertexAllocator>::apply(m_vertices, points);
        reserve_if_vector<face_type, FaceContainer, FaceAllocator>::apply(m_faces, 2 * points - 5);
    }

    triangulation(std::size_t points, std::size_t faces)
    {
        reserve_if_vector<vertex_type, VertexContainer, VertexAllocator>::apply(m_vertices, points);
        reserve_if_vector<face_type, FaceContainer, FaceAllocator>::apply(m_faces, faces);
    }

    template <typename InputIt>
    triangulation(InputIt begin, InputIt end)
    {
        m_vertices.assign(begin, end);
        reserve_if_vector<face_type, FaceContainer, FaceAllocator>
            ::apply(m_faces, 2 * m_vertices.size() - 5);
    }

    typename vertex_container::iterator vertices_begin()
    {
        return m_vertices.begin();
    }

    typename vertex_container::const_iterator vertices_begin() const
    {
        return m_vertices.cbegin();
    }

    typename vertex_container::iterator vertices_end()
    {
        return m_vertices.end();
    }

    typename vertex_container::const_iterator vertices_end() const
    {
        return m_vertices.cend();
    }

    vertex_iterator add_vertex(const Point& p)
    {
        return m_vertices.insert(m_vertices.end(), vertex_type{p, invalid});
    }

    face_container const& face_range() const
    {
        return m_faces;
    }

    vertex_container const& vertex_range() const
    {
        return m_vertices;
    }

    typename face_container::const_iterator faces_cbegin() const
    {
        return m_faces.cbegin();
    }

    typename face_container::iterator faces_begin()
    {
        return m_faces.begin();
    }

    typename face_container::const_iterator faces_cend() const
    {
        return m_faces.cend();
    }

    typename face_container::iterator faces_end()
    {
        return m_faces.end();
    }

    template <typename InputIt>
    void assign_vertices(InputIt begin, InputIt end)
    {
        m_vertices.assign(begin, end);
        m_faces.reserve(2 * m_vertices.size() - 5);
    }

    Point const& face_vertex(face_iterator f, face_vertex_index v)
    {
        return f -> m_v[v] -> m_p;
    }

    segment_type face_segment(halfedge_index e)
    {
        return segment_type(
            face_vertex(e.m_f, (e.m_v == 2 ? 0 : e.m_v +1)),
            face_vertex(e.m_f, (e.m_v == 0 ? 2 : e.m_v - 1)) );
    }

    Point& vertex(vertex_iterator v)
    {
        return v -> m_p;
    }

    face_iterator neighbour(face_iterator f, unsigned short v) const
    {
        return f -> m_f[v];
    }

    const_face_iterator neighbour(const_face_iterator f, unsigned short v) const
    {
        return f -> m_f[v];
    }

    face_vertex_index opposite(face_iterator f, unsigned short v) const
    {
        return f -> m_o[v];
    }

    halfedge_index opposite(halfedge_index const& e) const
    {
        return halfedge_index{ e.m_f -> m_f[e.m_v], e.m_f -> m_o[e.m_v] };
    }

    halfedge_index next(halfedge_index const& e) const
    {
        return halfedge_index{ e.m_f, static_cast<unsigned short>(e.m_v == 2 ? 0 : e.m_v + 1) };
    }

    halfedge_index prev(halfedge_index const& e) const
    {
        return halfedge_index{ e.m_f, static_cast<unsigned short>(e.m_v == 0 ? 2 : e.m_v - 1)};
    }

    vertex_iterator boundary_vertex() const
    {
        return m_boundary_vertex;
    }

    vertex_iterator boundary_next(vertex_iterator v) const
    {
        face_iterator fi = v -> m_f;
        if(fi -> m_v[0] == v )
            return fi -> m_v[1];
        else if(fi -> m_v[1] == v)
            return fi -> m_v[2];
        else
            return fi -> m_v[0];
    }

    vertex_iterator boundary_prev(vertex_iterator v) const
    {
        face_iterator fi = v -> m_f;
        unsigned short vi;
        if(fi -> m_v[0] == v) vi = 0;
        else if(fi -> m_v[1] == v) vi = 1;
        else vi = 2;
        if(m_faces.size()==1)
            return vi == 0 ? fi -> m_v[2] : fi -> m_v[vi-1];
        halfedge_index e = next(halfedge_index{fi, vi});
        while( opposite(e).m_f != invalid )
        {
            e = prev(opposite(e));
        }
        return e.m_f -> m_v[ e.m_v == 2 ? 0 : e.m_v + 1 ];
    }

    void clear()
    {
        m_vertices.clear();
        m_faces.clear();
    }

    std::size_t vertices() const
    {
        return m_vertices.size();
    }

    std::size_t faces() const
    {
        return m_faces.size();
    }

    void flip(const halfedge_index& e)
    {
        face_iterator fi1 = e.m_f;
        face_type& f1 = *fi1;
        face_iterator fi2 = f1.m_f[e.m_v];
        face_type& f2 = *fi2;
        unsigned short const& v1 = e.m_v;
        unsigned short const v2 = f1.m_o[v1];

        if( f1.m_v[ v1 == 0 ? 2 : v1 - 1 ] -> m_f == fi1 )
            f1.m_v[ v1 == 0 ? 2 : v1 - 1 ] -> m_f = fi2;
        if( f2.m_v[ v2 == 0 ? 2 : v2 - 1  ]-> m_f == fi2)
            f2.m_v[ v2 == 0 ? 2 : v2 - 1  ]-> m_f = fi1;
        f1.m_v[ v1 == 0 ? 2 : v1 - 1 ] = f2.m_v[ v2 ];
        f2.m_v[ v2 == 0 ? 2 : v2 - 1 ] = f1.m_v[ v1 ];

        f1.m_f[v1] = f2.m_f[ v2 == 2 ? 0 : v2 + 1 ];
        f1.m_o[v1] = f2.m_o[ v2 == 2 ? 0 : v2 + 1 ];
        if(f1.m_f[v1] != invalid) {
            f1.m_f[v1] -> m_f[f2.m_o[v2 == 2 ? 0 : v2 + 1]] = fi1;
            f1.m_f[v1] -> m_o[f2.m_o[v2 == 2 ? 0 : v2 + 1]] = v1;
        }
        f2.m_f[v2] = f1.m_f[ v1 == 2 ? 0 : v1 + 1 ];
        f2.m_o[v2] = f1.m_o[ v1 == 2 ? 0 : v1 + 1 ];
        if(f2.m_f[v2] != invalid) {
            f2.m_f[v2] -> m_f[f1.m_o[v1 == 2 ? 0 : v1 + 1]] = fi2;
            f2.m_f[v2] -> m_o[f1.m_o[v1 == 2 ? 0 : v1 + 1]] = v2;
        }
        f1.m_f[ v1 == 2 ? 0 : v1 + 1 ] = fi2;
        f1.m_o[ v1 == 2 ? 0 : v1 + 1 ] = v2 == 2 ? 0 : v2 + 1;
        f2.m_f[ v2 == 2 ? 0 : v2 + 1 ] = fi1;
        f2.m_o[ v2 == 2 ? 0 : v2 + 1 ] = v1 == 2 ? 0 : v1 + 1;
    }

    face_iterator add_face_on_boundary(halfedge_index e, vertex_iterator v)
    {
        const face_iterator f = e.m_f;
        const face_vertex_index adj = e.m_v;
        f -> m_o[adj] = 0;
        m_boundary_vertex = v;
        face_iterator pos = m_faces.insert(m_faces.end(),
            face_type{ {{v, f->m_v[ adj == 0 ? 2 : adj - 1 ], f->m_v[ adj == 2 ? 0 : adj + 1 ] }},
                {{f, invalid, invalid}},
                {{adj, 4, 4}}});
        f -> m_f[adj] = pos;
        v -> m_f = pos;
        m_faces.back().m_v[2]->m_f = pos;
        return pos;
    }

    fulledge_index next_around_vertex(fulledge_index e)
    {
        if(e.m_f2 == invalid) {
            face_vertex_index left_vi = (e.m_v1 == 0 ? 2 : e.m_v1 - 1);
            vertex_iterator left_v = e.m_f1->m_v[left_vi];
            face_iterator next_f = left_v->m_f;
            face_vertex_index next_vi = 
                (next_f->m_v[0] == left_v) ? 2 :
                ((next_f->m_v[1] == left_v) ? 0 : 1);
            return fulledge_index(invalid, 4, next_f, next_vi);
        } else {
            return fulledge_index(halfedge_index(e.m_f2, e.m_v2 == 0 ? 2 : e.m_v2 - 1));
        }
    }

    fulledge_index begin_vertex_edge(vertex_iterator vi)
    {
        face_iterator first_f = vi -> m_f;
        face_vertex_index next_vi = 
            (first_f->m_v[0] == vi) ? 2 :
            ((first_f->m_v[1] == vi) ? 0 : 1);
        return fulledge_index(invalid, 4, first_f, next_vi);
    }

    face_iterator add_isolated_face(vertex_iterator v1, vertex_iterator v2, vertex_iterator v3)
    {
        m_boundary_vertex = v1;
        face_iterator pos = m_faces.insert( m_faces.end(),
            face_type{
                {{ v1, v2, v3 }},
                {{ invalid, invalid, invalid }},
                {{4, 4, 4}} } );
        v1 -> m_f = v2 -> m_f = v3 -> m_f = pos;
        return pos;
    }

    halfedge_index face_edge(face_iterator f, face_vertex_index v = 0) const
    {
        return halfedge_index{f, v};
    }

    void connect(halfedge_index e1, halfedge_index e2)
    {
        face_iterator f1 = e1.m_f;
        face_iterator f2 = e2.m_f;
        unsigned short& v1 = e1.m_v;
        unsigned short& v2 = e2.m_v;
        f1 -> m_f[v1] = f2;
        f1 -> m_o[v1] = v2;
        f2 -> m_f[v2] = f1;
        f2 -> m_o[v2] = v1;
        if(f1 -> m_v[ (v1 == 2 ? 0 : v1 + 1) ] -> m_f == f1) {
            f1 -> m_v[ (v1 == 2 ? 0 : v1 + 1) ] -> m_f = f2;
        }
        if(f2 -> m_v[ (v2 == 2 ? 0 : v2 + 1) ] -> m_f == f2) {
            f2 -> m_v[ (v2 == 2 ? 0 : v2 + 1) ] -> m_f = f1;
        }
    }

/*
    void remove_face(face_index fi)
    {
        auto& f = m_faces[fi];
        if(m_faces.size() > 1)
            while(m_vertices[m_boundary_vertex].m_f == fi)
                m_boundary_vertex = boundary_next(m_boundary_vertex);
        for(face_vertex_index i = 0; i < 3 ; ++i)
            if(m_vertices[f.m_v[i]].m_f == fi)
                if(f.m_f[i == 0 ? 2 : i-1 ] != invalid)
                    m_vertices[f.m_v[i]].m_f = f.m_f[i == 0 ? 2 : i-1 ];
                else if (f.m_f[i == 2 ? 0 : i+1 ] != invalid)
                    m_vertices[f.m_v[i]].m_f = f.m_f[i == 2 ? 0 : i+1 ];
                else
                    m_vertices[f.m_v[i]].m_f = invalid;
        for(face_vertex_index i=0; i<3; ++i)
            m_faces[f.m_f[i]].m_f[f.m_o[i]] = invalid;
        for(face_vertex_index i=0; i<3; ++i)
            m_faces[f.m_f[i]].m_o[f.m_o[i]] = 4;
        for(face_index i = 0; i < m_faces.size(); ++i)
            for(face_vertex_index j = 0; j<3; ++j)
                if(m_faces[i].m_f[j] > fi && m_faces[i].m_f[j] != invalid)
                    m_faces[i].m_f[j]--;
        for(vertex_index i = 0; i < m_vertices.size(); ++i)
            if(m_vertices[i].m_f > fi && m_vertices[i].m_f != invalid)
                m_vertices[i].m_f--;
        m_faces.erase(m_faces.begin()+fi);
    }*/

    bool valid() const
    {
        bool valid = true;
        for(const_face_iterator fi = m_faces.begin(); fi != m_faces.end(); ++fi)
        {
            face_type const& f = *fi;
            for(unsigned short v = 0 ; v < 3 ; ++v)
            {
                if(f.m_f[v] == invalid)
                    continue;
                face_vertex_index const& o = f.m_o[v];
                valid = valid && (f.m_f[v] -> m_o[o] == v);
                if(!valid) {
                    std::cout << "1\n";
                    return false;
                }
                valid = valid && (f.m_f[v] -> m_f[o] == fi);
                if(!valid) {
                    std::cout << "2\n";
                    return false;
                }
                valid = valid && 
                    (f.m_v[ (v+1)%3 ] == f.m_f[v] -> m_v[ (o+2)%3 ]) 
                    && (f.m_v[ (v+2)%3 ] == f.m_f[v] -> m_v[ (o+1)%3 ]);
                if(!valid) {
                    std::cout << "3\n";
                    return false;
                }
                if(f.m_o[v] == 4) {
                    unsigned short next = (v + 1) % 3;
                    valid == valid && f.m_v[next]->m_f == fi;
                    if(!valid) {
                        std::cout << "6\n";
                        return false;
                    }

                }
            }
            valid = valid && 
                (strategy::side::side_by_triangle<>
                    ::apply(f.m_v[0]->m_p, f.m_v[1]->m_p, f.m_v[2]->m_p) > 0);
            if(!valid) {
                std::cout << "4\n";
                return false;
            }
        }
        for(const_vertex_iterator vi = m_vertices.cbegin(); vi != m_vertices.cend(); ++vi)
        {
            vertex_type const& v = *vi;
            if(v.m_f == invalid) continue;
            bool found = false;
            for(face_vertex_index vj = 0; vj < 3 ; ++vj)
            {
                found = found || (v.m_f -> m_v[vj] == vi);
            }
            valid = valid && found;
            if(!valid) {
                std::cout << "Error 5 at " 
                    << std::distance<const_vertex_iterator>(m_vertices.cbegin(), vi) << "\n";
                return false;
            }
        }
        return valid;
    }
private:
    vertex_container m_vertices;
    face_container m_faces;
    vertex_iterator m_boundary_vertex;
};

template< typename Point >
struct edge_ref
{
    typename triangulation<Point>::halfedge_index m_e;
    triangulation<Point>& m_t;
};

template< typename Point >
using triangulation_face_range = typename triangulation<Point>::face_container;

template< typename Point >
using triangulation_vertex_range = typename triangulation<Point>::vertex_container;

} // namespace model

/*
template<typename Point>
struct triangulation_face<model::triangulation<Point>>
{
    typename model::triangulation<Point>::face_type const& 
        get(model::triangulation<Point> const& t, typename model::triangulation<Point>::face_index const& fi)
    {
        return t.face(fi);
    }
};*/

struct triangulation_tag {};

#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{

template <typename Point>
struct tag< model::triangulation<Point> >
{
    typedef triangulation_tag type;
};

template <typename Point>
struct point_type< model::triangulation<Point> >
{
    typedef typename model::triangulation<Point>::point_type type;
};

template<typename Triangulation>
struct dimension< model::vertex_ref<Triangulation> > : boost::mpl::int_<2> {};

template<typename Triangulation>
struct tag< model::vertex_ref<Triangulation> >
{ typedef point_tag type; };

template<typename Triangulation>
struct coordinate_type< model::vertex_ref<Triangulation> >
{ typedef typename coordinate_type<typename Triangulation::point_type>::type type; };

template<typename Triangulation>
struct coordinate_system< model::vertex_ref<Triangulation> >
{ typedef typename coordinate_system<typename Triangulation::point_type>::type type; };

template<typename Triangulation, std::size_t Dimension>
struct access<model::vertex_ref<Triangulation>, Dimension>
{
    static typename coordinate_type<typename Triangulation::point_type>::type
        get(model::vertex_ref<Triangulation> const& p)
    {
        return boost::geometry::get<Dimension>(p.m_p);
    }
};

template<typename Point> struct tag< model::edge_ref<Point> >
{ typedef segment_tag type; };

template<typename Point> struct point_type< model::edge_ref<Point> >
{ typedef Point type; };

template<typename Point, std::size_t Dimension>
struct indexed_access<model::edge_ref<Point>, 0, Dimension>
{
    static typename coordinate_type<Point>::type get(model::edge_ref<Point> const& p)
    {
        return get<Dimension>(p.m_t.face_vertex(p.m_e.m_f, (p.m_e.m_v == 2 ? 0 : p.m_e.m_v + 1)));
    }
};

template<typename Point, std::size_t Dimension>
struct indexed_access<model::edge_ref<Point>, 1, Dimension>
{
    static typename coordinate_type<Point>::type get(model::edge_ref<Point> const& p)
    {
        return get<Dimension>(p.m_t.face_vertex(p.m_e.m_f, (p.m_e.m_v == 0 ? 2 : p.m_e.m_v - 1)));
    }
};

template<typename Triangulation> struct tag<model::face_ref<Triangulation>>
{ typedef ring_tag type; };

template<typename Triangulation> struct point_order<model::face_ref<Triangulation>>
{ static const order_selector value = counterclockwise; };

template<typename Triangulation>
struct closure<model::face_ref<Triangulation>>
{
    static const closure_selector value = open;
};

} // namespace traits
#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS

template<typename Point>
struct face_range_type<model::triangulation<Point>> {
    typedef typename model::triangulation_face_range<Point> type;
};

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
inline typename model::triangulation
<
    Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator
>::vertex_container const&
    vertex_range(
        model::triangulation<Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
            const& t)
{ return t.vertex_range(); }

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
inline typename model::triangulation
<
    Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator
>
    ::face_container const&
    face_range(model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> const& t)
{ return t.face_range(); }

template
<
    typename Iterator
>
struct indirect_range
{
    std::vector<Iterator> _iterators;
public:
    void push_back(Iterator i) {_iterators.push_back(i);}
    typedef boost::indirect_iterator<typename std::vector<Iterator>::const_iterator> const_iterator;
    typedef boost::indirect_iterator<typename std::vector<Iterator>::const_iterator> iterator;

    const_iterator begin() const { return const_iterator(_iterators.cbegin()); }
    const_iterator end() const { return const_iterator(_iterators.cend()); }
    iterator begin() { return iterator(_iterators.begin()); }
    iterator end() { return iterator(_iterators.end()); }
};

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
inline indirect_range<
    typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
        ::face_iterator>
    face_adjacent_range(
        model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>& t,
        typename model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
            ::const_face_iterator fi
    )
{
    typedef typename model::triangulation
    <
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator
    >
        triangulation;
    typedef typename triangulation::face_iterator face_iterator;
    typedef typename triangulation::face_type face_type;
    face_type const& f = *fi;
    face_iterator const invalid = t.invalid;
    indirect_range<face_iterator> out;
    if(f.m_f[0] != invalid ) out.push_back(f.m_f[0]);
    if(f.m_f[1] != invalid ) out.push_back(f.m_f[1]);
    if(f.m_f[2] != invalid ) out.push_back(f.m_f[2]);
    return out;
}

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator,
    typename Iterator
>
inline indirect_range<
    typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>::face_iterator>
    face_adjacent_range(
        model::triangulation<Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> & t,
        boost::indirect_iterator<Iterator> fi
    )
{
    return face_adjacent_range(t, *fi.base());
}

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
inline indirect_range<typename model::triangulation<
    Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>::face_iterator
>
    face_incident_range(
        model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> & t,
        typename model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
            ::face_iterator fi
    )
{
    typedef typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> 
        triangulation;
    typedef typename triangulation::face_iterator face_iterator;
    typedef typename triangulation::vertex_iterator vertex_iterator;
    typedef typename triangulation::face_vertex_index face_vertex_index;
    typedef typename triangulation::face_type face_type;
    typedef typename triangulation::const_face_iterator const_face_iterator;
    face_type const& f = *fi;
    indirect_range<face_iterator> out;
    face_iterator const invalid = t.invalid;
    for(face_vertex_index i = 0; i < 3; ++i)
    {
        face_iterator n = t.neighbour(fi, i);
        face_iterator m = t.neighbour(fi, (i == 2 ? 0 : i + 1));
        face_iterator f_prev = fi;
        face_vertex_index v_prev = i;
        if(n != invalid) {
            f_prev = n;
            v_prev = fi -> m_o[i];
            v_prev = (v_prev == 0 ? 2 : v_prev - 1);
        }
        while(true)
        {
            face_iterator next = t.neighbour(f_prev, v_prev);
            if(next == invalid) {
                face_vertex_index j = (i == 2 ? 0 : i + 1);
                face_iterator m = t.neighbour(fi, j);
                if(m == invalid) break;
                face_vertex_index prev_vertex_index = (i == 0 ? 2 : i - 1 );
                vertex_iterator const& prev_vertex_it = f.m_v[prev_vertex_index];
                face_type& first = *prev_vertex_it->m_f;
                out.push_back(prev_vertex_it->m_f);
                f_prev = next = prev_vertex_it->m_f;
                if(first.m_v[0] == prev_vertex_it) v_prev = 1;
                else if(first.m_v[1] == prev_vertex_it) v_prev = 2;
                else v_prev = 0;
                continue;
            } else {
                if(next == fi) break;
                out.push_back(next);
                v_prev = f_prev -> m_o[v_prev];
                v_prev = (v_prev == 0 ? 2 : v_prev - 1);
                f_prev = next;
            }
        }
    }
    return out;
}

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator,
    typename Iterator
>
inline indirect_range<
    typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
        ::face_iterator>
    face_incident_range(
        model::triangulation<Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> & t,
        boost::indirect_iterator<Iterator> fi
    )
{
    return face_incident_range(t, *fi.base());
}

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
inline indirect_range<
    typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
        ::face_iterator>
    face_incident_range(
        model::triangulation<Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> & t,
        typename model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
            ::vertex_iterator vi
    )
{
    typedef typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> 
        triangulation;
    typedef typename triangulation::face_iterator face_iterator;
    typedef typename triangulation::fulledge_index fulledge_index;
    indirect_range<face_iterator> out;
    fulledge_index e = t.begin_vertex_edge(vi);
    face_iterator first_face = e.m_f2;
    out.push_back(first_face);
    while(true) {
        e = t.next_around_vertex(e);
        if(e.m_f2 == t.invalid || e.m_f2 == first_face) break;
        out.push_back(e.m_f2);
    }
    return out;
}

template
<
    typename Point,
    template<typename, typename> class VertexContainer,
    template<typename, typename> class FaceContainer,
    template<typename> class VertexAllocator,
    template<typename> class FaceAllocator
>
inline indirect_range<typename model::triangulation
    <Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>::vertex_iterator>
    vertex_incident_range(
        model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator> & t,
        typename model::triangulation<
            Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
            ::vertex_iterator vi
    )
{
    typedef typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
        ::face_iterator face_iterator;
    typedef typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
        ::vertex_iterator vertex_iterator;
    typedef typename model::triangulation<
        Point, VertexContainer, FaceContainer, VertexAllocator, FaceAllocator>
        ::fulledge_index fulledge_index;
    indirect_range<vertex_iterator> out;
    fulledge_index e = t.begin_vertex_edge(vi);
    face_iterator first_face = e.m_f2;
    out.push_back(e.m_f2 -> m_v[ e.m_v2 == 0 ? 2 : e.m_v2 - 1 ]);
    while(true) {
        e = t.next_around_vertex(e);
        if(e.m_f2 == t.invalid || e.m_f2 == first_face) break;
        out.push_back(e.m_f2 -> m_v[ e.m_v2 == 0 ? 2 : e.m_v2 - 1 ]);
    }
    return out;
}

} // namespace geometry

} // namespace boost

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_GEOMETRIES_TRIANGULATION_HPP
