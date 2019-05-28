#ifndef BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_GEOMETRIES_TRIANGULATION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_GEOMETRIES_TRIANGULATION_HPP

#include <vector>

#include <boost/range.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>

#include <boost/geometry/strategies/cartesian/side_by_triangle.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/views/detail/points_view.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>

#include <boost/geometry/extensions/triangle/triangle.hpp>

namespace boost { namespace geometry
{

namespace model
{

template <typename Point>
class triangulation;

template<typename Point>
struct face_ref
{
    std::array<std::size_t, 3> m_v;
    std::array<std::size_t, 3> m_f;
    std::array<unsigned short, 3> m_o;
    const triangulation<Point>* _t;

private:
    struct face_vertex_iterator : public boost::iterator_facade
        <
            face_vertex_iterator,
            Point const,
            boost::random_access_traversal_tag
        >
    {
       //Constructor: Begin iterator
        inline face_vertex_iterator(face_ref const* f) : m_f(f),m_v(0) {}
        //Constructor: End iterator
        inline face_vertex_iterator(face_ref const* f, bool) : m_f(f), m_v(3) {}
        inline face_vertex_iterator() : m_f(nullptr), m_v(3) {}

        typedef short difference_type;
    private:
        friend class boost::iterator_core_access;

        inline Point const& dereference() const
        {
            return (m_f->_t->vertices_begin() + m_f->m_v[m_v])->m_p;
        }

        inline bool equal(face_vertex_iterator const& other) const
        {
            return other.m_f == this->m_f && other.m_v == this->m_v;
        }

        inline void increment()
        {
            m_v++;
        }

        inline void decrement()
        {
            m_v--;
        }

        inline difference_type distance_to(face_vertex_iterator const& other) const
        {
           return other.m_v - this->m_v;
        }

        inline void advance(difference_type n)
        {
            m_v += n;
        }
        face_ref const* m_f;
        short m_v;
    };
public:
    typedef face_vertex_iterator const_iterator;
    typedef face_vertex_iterator iterator; // must be defined

    const_iterator begin() const { return const_iterator(this); }
    const_iterator end() const { return const_iterator(this, true); }
};

template<typename Point>
struct vertex_ref
{
    Point m_p;
    std::size_t m_f;
};

template <typename Point>
class triangulation
{
private:
    BOOST_CONCEPT_ASSERT( (concepts::Point<Point>) );
public:

    typedef std::size_t face_index;
    typedef std::size_t vertex_index;
    typedef unsigned short face_vertex_index;
	struct edge_index
	{
        edge_index(face_index f, face_vertex_index v):m_f(f), m_v(v) {}
		face_index m_f;
		face_vertex_index m_v;
	};

    void debug_print() 
    {
        std::cout << "Vertices: \n";
        for(std::size_t i = 0; i < m_vertices.size(); ++i)
            std::cout << "Vertex " << i << ": ( " << get<0>(m_vertices[i].m_p) << " , " << get<1>(m_vertices[i].m_p) << "), touches face: " << m_vertices[i].m_f << "\n";
        std::cout << "Faces: \n";
        for(std::size_t i = 0; i < m_faces.size(); ++i) {
            debug_print_face(i);
        }
        std::cout << "boundary vertex: " << m_boundary_vertex << "\n";
    }

    void debug_print_face(std::size_t i) {
            const auto& f = m_faces[i];
            std::cout << "Face " << i << ":\n";
            for(unsigned short v = 0; v < 3; ++v)
                std::cout << "Vertex " << v << ": " << f.m_v[v] << ", Neighbour: " << f.m_f[v] << ", Opposite: " << f.m_o[v] << "\n";
    }

	typedef typename std::vector<vertex_ref<Point>> vertex_container;
	typedef typename std::vector<face_ref<Point>> face_container;
    typedef typename coordinate_type<Point>::type coordinate_type;
    typedef typename model::segment<Point> segment_type;
    typedef Point point_type;
	static const face_index invalid = -1;
    triangulation() {}
    triangulation(std::size_t points)
	{
		m_vertices.reserve(points);
		m_faces.reserve(2 * points - 5);
	}

    triangulation(std::size_t points, std::size_t faces)
	{
		m_vertices.reserve(points);
		m_faces.reserve(faces);
	}

	template <typename InputIt>
	triangulation(InputIt begin, InputIt end)
	{
		m_vertices.assign(begin, end);
		m_faces.reserve(2 * m_vertices.size() - 5);
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

	vertex_index add_vertex(const Point& p)
	{
        m_vertices.push_back(vertex_ref<Point>{p,std::size_t(-1)});
		return m_vertices.size()-1;
	}

    face_container const& face_range() const
    {
        return m_faces;
    }

    vertex_container const& vertex_range() const
    {
        return m_vertices;
    }

    typename face_container::const_iterator faces_begin() const
	{
		return m_faces.cbegin();
	}

    typename face_container::const_iterator faces_end() const
	{
		return m_faces.cend();
	}

    template <typename InputIt>
    void assign_vertices(InputIt begin, InputIt end)
    {
        m_vertices.assign(begin, end);
        m_faces.reserve(2 * m_vertices.size() - 5);
    }

    Point& face_vertex(face_index f, face_vertex_index v)
    {
        return m_vertices[m_faces[f].m_v[v]].m_p;
    }

    segment_type face_segment(edge_index e)
    {
        return segment_type(face_vertex(e.m_f, (e.m_v == 2 ? 0 : e.m_v +1)), face_vertex(e.m_f, (e.m_v == 0 ? 2 : e.m_v - 1)) );
    }

    Point& vertex(vertex_index v)
    {
        return m_vertices[v].m_p;
    }

    face_index neighbour(face_index f, unsigned short v)
    {
        return m_faces[f].m_f[v];
    }

    face_vertex_index opposite(face_index f, unsigned short v)
    {
        return m_faces[f].m_o[v];
    }

    edge_index opposite(edge_index const& e)
    {
        return edge_index{ m_faces[e.m_f].m_f[e.m_v], m_faces[e.m_f].m_o[e.m_v] };
    }

    edge_index next(edge_index const& e)
    {
        return edge_index{ e.m_f, static_cast<unsigned short>(e.m_v == 2 ? 0 : e.m_v + 1) };
    }

    edge_index prev(edge_index const& e)
    {
        return edge_index{ e.m_f, static_cast<unsigned short>(e.m_v == 0 ? 2 : e.m_v - 1)};
    }

    vertex_index boundary_vertex()
    {
        return m_boundary_vertex;
    }

    vertex_index boundary_next(vertex_index const& v)
    {
        face_index fi = m_vertices[v].m_f;
        unsigned short vi;
        if(m_faces[fi].m_v[0] == v) vi = 0;
        else if(m_faces[fi].m_v[1] == v) vi = 1;
        else vi = 2;
        if(m_faces.size()==1) return vi == 2 ? m_faces[fi].m_v[0] : m_faces[fi].m_v[vi+1];
        edge_index e = prev(edge_index{fi, vi});
        while( opposite(e).m_f != invalid )
        {
            e = next(opposite(e));
        }
        return m_faces[e.m_f].m_v[ e.m_v == 0 ? 2 : e.m_v - 1 ];
    }

    vertex_index boundary_prev(vertex_index const& v)
    {
        face_index fi = m_vertices[v].m_f;
        unsigned short vi;
        if(m_faces[fi].m_v[0] == v) vi = 0;
        else if(m_faces[fi].m_v[1] == v) vi = 1;
        else vi = 2;
        if(m_faces.size()==1) return vi == 0 ? m_faces[fi].m_v[2] : m_faces[fi].m_v[vi-1];
        edge_index e = next(edge_index{fi, vi});
        while( opposite(e).m_f != invalid )
        {
            e = prev(opposite(e));
        }
        return m_faces[e.m_f].m_v[ e.m_v == 2 ? 0 : e.m_v + 1 ];
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
	
    void flip(const edge_index& e)
	{
        std::size_t const& fi1 = e.m_f;
		face_ref<Point>& f1 = m_faces[fi1];
        std::size_t const fi2 = f1.m_f[e.m_v];
		face_ref<Point>& f2 = m_faces[fi2];
        unsigned short const& v1 = e.m_v;
        unsigned short const v2 = f1.m_o[v1];

        m_vertices[ f1.m_v[ v1 == 0 ? 2 : v1 - 1 ] ].m_f = fi2;
        m_vertices[ f2.m_v[ v2 == 0 ? 2 : v2 - 1 ] ].m_f = fi1;
		f1.m_v[ v1 == 0 ? 2 : v1 - 1 ] = f2.m_v[ v2 ];
		f2.m_v[ v2 == 0 ? 2 : v2 - 1 ] = f1.m_v[ v1 ];

		f1.m_f[v1] = f2.m_f[ v2 == 2 ? 0 : v2 + 1 ];
        f1.m_o[v1] = f2.m_o[ v2 == 2 ? 0 : v2 + 1 ];
        if(f1.m_f[v1] != invalid) {
            m_faces[f1.m_f[v1]].m_f[f2.m_o[v2 == 2 ? 0 : v2 + 1]] = fi1;
            m_faces[f1.m_f[v1]].m_o[f2.m_o[v2 == 2 ? 0 : v2 + 1]] = v1;
        }
        f2.m_f[v2] = f1.m_f[ v1 == 2 ? 0 : v1 + 1 ];
        f2.m_o[v2] = f1.m_o[ v1 == 2 ? 0 : v1 + 1 ];
        if(f2.m_f[v2] != invalid) {
            m_faces[f2.m_f[v2]].m_f[f1.m_o[v1 == 2 ? 0 : v1 + 1]] = fi2;
            m_faces[f2.m_f[v2]].m_o[f1.m_o[v1 == 2 ? 0 : v1 + 1]] = v2;
        }
        f1.m_f[ v1 == 2 ? 0 : v1 + 1 ] = fi2;
        f1.m_o[ v1 == 2 ? 0 : v1 + 1 ] = v2 == 2 ? 0 : v2 + 1;
		f2.m_f[ v2 == 2 ? 0 : v2 + 1 ] = fi1;
        f2.m_o[ v2 == 2 ? 0 : v2 + 1 ] = v1 == 2 ? 0 : v1 + 1;
	}

	face_index add_face_on_boundary(edge_index e, vertex_index v)
	{
        const auto f = e.m_f;
        const auto adj = e.m_v;
		m_faces[f].m_f[adj] = m_faces.size();
        m_faces[f].m_o[adj] = 0;
        m_vertices[v].m_f = m_faces.size();
        m_boundary_vertex = v;
		m_faces.push_back(face_ref<Point>{ {{v, m_faces[f].m_v[ adj == 0 ? 2 : adj - 1 ], m_faces[f].m_v[ adj == 2 ? 0 : adj + 1 ] }},
			{{f, invalid, invalid}},
            {{adj, 4, 4}},
            this});
		return m_faces.size()-1;
	}

	face_index add_isolated_face(vertex_index v1, vertex_index v2, vertex_index v3)
	{
        m_boundary_vertex = v1;
        m_vertices[v1].m_f = m_vertices[v2].m_f = m_vertices[v3].m_f = m_faces.size();
        m_faces.insert( m_faces.end(), face_ref<Point>{ {{ v1, v2, v3 }}, {{ invalid, invalid, invalid }}, {{4, 4, 4}}, this } );
		return m_faces.size() - 1;
	}

    edge_index face_edge(face_index f, face_vertex_index v = 0)
    {
        return edge_index{f, v};
    }

    void connect(edge_index e1, edge_index e2)
    {
        face_index const& f1 = e1.m_f;
        face_index const& f2 = e2.m_f;
        unsigned short& v1 = e1.m_v;
        unsigned short& v2 = e2.m_v;
        m_faces[f1].m_f[v1] = f2;
        m_faces[f1].m_o[v1] = v2;
        m_faces[f2].m_f[v2] = f1;
        m_faces[f2].m_o[v2] = v1;
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

    bool valid()
    {
        bool valid = true;
        for(std::size_t fi = 0; fi < m_faces.size(); ++fi)
        {
            auto const& f = m_faces[fi];
            for(unsigned short v = 0 ; v < 3 ; ++v)
            {
                if(f.m_f[v] == invalid)
                    continue;
                auto const& o = f.m_o[v];
                valid = valid && (m_faces[f.m_f[v]].m_o[o] == v);
                if(!valid) {
                    std::cout << "1\n";
                    return false;
                }
                valid = valid && (m_faces[f.m_f[v]].m_f[o] == fi);
                if(!valid) {
                    std::cout << "2\n";
                    return false;
                }
                valid = valid && (f.m_v[ (v+1)%3 ] == m_faces[f.m_f[v]].m_v[ (o+2)%3 ]) && (f.m_v[ (v+2)%3 ] == m_faces[f.m_f[v]].m_v[ (o+1)%3 ]);
                if(!valid) {
                    std::cout << "3\n";
                    return false;
                }
            }
            valid = valid && (strategy::side::side_by_triangle<>::apply(m_vertices[f.m_v[0]].m_p, m_vertices[f.m_v[1]].m_p, m_vertices[f.m_v[2]].m_p) > 0);
            if(!valid) {
                std::cout << "4\n";
                return false;
            }
        }
        for(std::size_t vi = 0; vi < m_vertices.size(); ++vi)
        {
            auto const& v = m_vertices[vi];
            if(v.m_f == invalid) continue;
            bool found = false;
            for(unsigned short vj = 0; vj < 3 ; ++vj)
            {
                found = found || (m_faces[v.m_f].m_v[vj] == vi);
            }
            valid = valid && found;
            if(!valid) {
                std::cout << "Error 5 at " << vi << "\n";
                return false;
            }
        }
        return valid;
    }
private:
    vertex_container m_vertices;
    face_container m_faces;
    std::size_t m_boundary_vertex;
};

template< typename Point >
struct edge_ref
{
    typename triangulation<Point>::edge_index m_e;    
    triangulation<Point>& m_t;
};

template< typename Point >
using triangulation_face_range = typename triangulation<Point>::face_container;

template< typename Point >
using triangulation_vertex_range = typename triangulation<Point>::vertex_container;

} // namespace model

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

template<typename Point> struct dimension< model::vertex_ref<Point> > : boost::mpl::int_<2> {};

template<typename Point> struct tag< model::vertex_ref<Point> >
{ typedef point_tag type; };

template<typename Point> struct coordinate_type< model::vertex_ref<Point> >
{ typedef typename coordinate_type<Point>::type type; };

template<typename Point> struct coordinate_system< model::vertex_ref<Point> >
{ typedef typename coordinate_system<Point>::type type; };

template<typename Point, std::size_t Dimension>
struct access<model::vertex_ref<Point>, Dimension>
{
    static typename coordinate_type<Point>::type get(model::vertex_ref<Point> const& p)
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
        return get<Dimension>( p.m_t.face_vertex(p.m_e.m_f, (p.m_e.m_v == 2 ? 0 : p.m_e.m_v + 1) ) );
    }                                                 
};

template<typename Point, std::size_t Dimension>
struct indexed_access<model::edge_ref<Point>, 1, Dimension>
{
    static typename coordinate_type<Point>::type get(model::edge_ref<Point> const& p)  
    {
        return get<Dimension>( p.m_t.face_vertex(p.m_e.m_f, (p.m_e.m_v == 0 ? 2 : p.m_e.m_v - 1) ) );  
    }                                                 
};

template<typename Point> struct tag<model::face_ref<Point>>
{ typedef ring_tag type; };

template<typename Point> struct point_order<model::face_ref<Point>>
{ static const order_selector value = counterclockwise; };

template<typename Point>
struct closure<model::face_ref<Point>>
{
    static const closure_selector value = open;
};

} // namespace traits
#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS

template<typename Triangulation>
struct face_range_type {};

template<typename Point>
struct face_range_type<model::triangulation<Point>> {
    typedef typename model::triangulation_face_range<Point> type;
};

template<typename Triangulation>
typename face_range_type<Triangulation>::type face_range(Triangulation const& t);

template<typename Point>
model::triangulation_vertex_range<Point> vertex_range(model::triangulation<Point> const& t)
{ return t.vertex_range(); }

template<typename Point>
model::triangulation_face_range<Point> face_range(model::triangulation<Point> const& t)
{ return t.face_range(); }

} // namespace geometry

} // namespace boost

#endif // BOOST_GEOMETRY_EXTENSIONS_TRIANGULATION_GEOMETRIES_TRIANGULATION_HPP
