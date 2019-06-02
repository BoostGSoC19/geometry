// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2019 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_RANDOM_DETAIL_UNIFORM_POINT_DISTRIBUTION_HPP
#define BOOST_GEOMETRY_EXTENSIONS_RANDOM_DETAIL_UNIFORM_POINT_DISTRIBUTION_HPP

#include <boost/range.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/assert.hpp>
#include <boost/geometry/util/for_each_coordinate.hpp>
#include <boost/geometry/util/select_most_precise.hpp>
#include <boost/geometry/views/box_view.hpp>
#include <boost/geometry/views/segment_view.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/algorithms/append.hpp>
#include <boost/geometry/algorithms/envelope.hpp>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>

namespace boost { namespace geometry { namespace random { namespace detail
{

    template<typename Point, typename DomainGeometry>
    class uniform_point_distribution_base
    {
    public:
        typedef Point result_type;
        class param_type {
        public:
            param_type(DomainGeometry const& domain):_domain(domain) {}
            param_type(param_type const& p):_domain(p._domain) {}
            DomainGeometry const& domain() const {return _domain;}
            bool operator==(param_type const& rhs) {return equals(_domain,rhs._domain);}
        private:
            DomainGeometry _domain;
        };
        uniform_point_distribution_base():_param(DomainGeometry()) {}
        uniform_point_distribution_base(param_type const& param):_param(param) {}
        uniform_point_distribution_base(DomainGeometry const& domain):_param(param_type(domain)) {}
        void reset() {}
        param_type param() const {return _param;}
        DomainGeometry const& domain() const {return _param.domain();}
        bool operator==(uniform_point_distribution_base const& rhs) {return _param==rhs._param;}
        void param(param_type const& p) {_param = p;}

    protected:
        param_type _param;
    };

    template
    <
        typename Point,
        typename DomainGeometry,
        typename DomainTag,
        typename SingleOrMulti,
        typename CSTag,
        std::size_t dimension
    >
    class uniform_point_distribution {};

    template<typename Point, typename DomainGeometry, typename CSTag, std::size_t dimension>
    class uniform_point_distribution<Point, DomainGeometry, pointlike_tag, single_tag, CSTag, dimension> 
        : public uniform_point_distribution_base<Point, DomainGeometry>
    {
        typedef uniform_point_distribution_base<Point, DomainGeometry> base;
    public:
        using base::base;
        template<typename Generator>
        Point operator()(Generator& gen, typename base::param_type const& p) {return this->domain();}
        template<typename Generator>
        Point operator()(Generator& gen) {return (*this)(gen, this->_param);}
    };

    template<typename Point, typename DomainGeometry, typename CSTag, std::size_t dimension>
    class uniform_point_distribution<Point, DomainGeometry, pointlike_tag, multi_tag, CSTag, dimension>
        : public uniform_point_distribution_base<Point, DomainGeometry>
    {
        typedef uniform_point_distribution_base<Point, DomainGeometry> base;
    public:
        using base::base;
        template<typename Generator>
        Point operator()(Generator& gen, typename base::param_type const& p) 
        {
            std::size_t count = num_points(p.domain());
            std::uniform_int_distribution<std::size_t> dis(0, count-1);
            return *(boost::begin(p.domain())+dis(gen));
        }
        template<typename Generator>
        Point operator()(Generator& gen) {return (*this)(gen, this->_param);}
    };

    template
    <
        typename Point,
        typename Box,
        typename Generator,
        bool isIntegral = boost::is_integral<typename coordinate_type<Point>::type>::type::value
    >
    struct interval_sample {};

    template<typename Point, typename Box, typename Generator>
    struct interval_sample<Point, Box, Generator, true>
    {
        Box const& b;
        Generator& gen;
        inline interval_sample(Box const& b, Generator& gen):b(b),gen(gen) {}

        template <typename PointDst, std::size_t Index>
        inline void apply(PointDst& point_dst) const
        {
            std::uniform_int_distribution<typename coordinate_type<Point>::type>
                dist(get<min_corner,Index>(b),get<max_corner,Index>(b));
            set<Index>(point_dst,dist(gen));
        }
    };

    template<typename Point, typename Box, typename Generator>
    struct interval_sample<Point, Box, Generator, false>
    {
        Box const& b;
        Generator& gen;
        inline interval_sample(Box const& b, Generator& gen):b(b),gen(gen) {}

        template <typename PointDst, std::size_t Index>
        inline void apply(PointDst& point_dst) const
        {
            std::uniform_real_distribution<typename coordinate_type<Point>::type>
                dist(get<min_corner,Index>(b),get<max_corner,Index>(b));
            set<Index>(point_dst,dist(gen));
        }
    };

    template<typename Point, typename DomainGeometry, std::size_t dimension>
    class uniform_point_distribution<Point, DomainGeometry, box_tag, single_tag, cartesian_tag, dimension> 
        : public uniform_point_distribution_base<Point, DomainGeometry>
    {
        typedef uniform_point_distribution_base<Point, DomainGeometry> base;
    public:
        using base::base;
        template<typename Generator>
        Point operator()(Generator& gen, typename base::param_type const& p) 
        {
            model::box<Point> cached_box;
            assign(cached_box, p.domain());
            Point out;
            for_each_coordinate(out, interval_sample<Point, DomainGeometry, Generator>(p.domain(), gen));
            return out;
        }
        template<typename Generator>
        Point operator()(Generator& gen)
        {
            return (*this)(gen, this->_param);
        }
    };

    template<typename Point, typename PointIn, typename LengthType>
    Point sample_segment(PointIn const& p1, PointIn const& p2, LengthType const& r)
    {
        Point out;
        assign(out,p2);
        subtract_point(out, p1);
        multiply_value(out, r);
        add_point(out, p1);
        return out;
    }

    template<typename Point, typename PointVec, typename IndexVec, typename LengthVec, typename LengthType>
    Point sample_multi_line(PointVec const& point_cache, IndexVec const& skip_list, LengthVec const& accumulated_lengths, LengthType const& r)
    {
        std::size_t i = std::distance(accumulated_lengths.begin(),
            std::lower_bound(accumulated_lengths.begin(), accumulated_lengths.end(), r));
        std::size_t offset = std::distance(skip_list.begin(),
            std::lower_bound(skip_list.begin(), skip_list.end(), i));
        return sample_segment<Point>(point_cache[i+offset-1], point_cache[i+offset],
            (r-accumulated_lengths[i-1])/(accumulated_lengths[i]-accumulated_lengths[i-1]));
    }

    template<typename Point, typename DomainGeometry, typename SingleOrMulti, typename CSTag, std::size_t dimension>
    class uniform_point_distribution<Point, DomainGeometry, areal_tag, SingleOrMulti, CSTag, dimension>
        : public uniform_point_distribution_base<Point, DomainGeometry>
    {
        typedef uniform_point_distribution_base<Point, DomainGeometry> base;
    public:
        typedef typename uniform_point_distribution_base<Point, DomainGeometry>::param_type param_type;
        using base::base;
        using base::param;
        uniform_point_distribution() {}
        uniform_point_distribution(DomainGeometry const& g):base(g) {envelope(g,_cached_box);}
        void param(param_type const& p) {this->_param=p; envelope(p.domain(),_cached_box);}
        template<typename Generator>
        Point operator()(Generator& gen, typename base::param_type const& p) 
        {
            model::box<Point> _cached_box;
            envelope(p.domain(), _cached_box);
            uniform_point_distribution<Point, model::box<Point>, box_tag, single_tag, CSTag, dimension>
                box_dist(_cached_box);
            Point out = box_dist(gen);
            while(!within(out, p.domain()))
                out = box_dist(gen);
            return out;
        }
        template<typename Generator>
        Point operator()(Generator& gen)
        {
            uniform_point_distribution<Point, model::box<Point>, box_tag, single_tag, CSTag, dimension>
                box_dist(_cached_box);
            Point out = box_dist(gen);
            while(!within(out, this->_param.domain()))
                out = box_dist(gen);
            return out;
        }
    private:
        model::box<Point> _cached_box;
    };

    template<typename Point, typename DomainGeometry, typename SingleOrMulti, std::size_t dimension>
    class uniform_point_distribution<Point, DomainGeometry, linear_tag, SingleOrMulti, cartesian_tag, dimension>
        : public uniform_point_distribution_base<Point, DomainGeometry>
    {
        typedef uniform_point_distribution_base<Point, DomainGeometry> base;
        typedef typename default_length_result<DomainGeometry>::type length_type;
        typedef typename point_type<DomainGeometry>::type domain_point_type;
    public:
        using base::base;
        using base::param;
        uniform_point_distribution(DomainGeometry const& g):base(g) {
            init(skip_list, point_cache, accumulated_lengths, this->_param);
        }
        void init(std::vector<std::size_t>& skip_list,
            std::vector<domain_point_type>& point_cache,
            std::vector<length_type>& accumulated_lengths,
            typename base::param_type const& p)
        {
            std::size_t i = 0;
            point_cache.push_back(*segments_begin(p.domain())->first);
            accumulated_lengths.push_back(0);
            for(auto it = segments_begin(p.domain()); it!=segments_end(p.domain()); ++it) {
                accumulated_lengths.push_back(accumulated_lengths.back()+length(*it));
                if(!equals(point_cache.back(),*it->first)) {
                    point_cache.push_back(*it->first);
                    skip_list.push_back(i);
                }
                point_cache.push_back(*it->second);
                ++i;
            }
        }
        void param(typename base::param_type const& p)
        {
            this->_param = p;
            skip_list.clear();
            point_cache.clear();
            accumulated_lengths.clear();
            init(skip_list, point_cache, accumulated_lengths, p);
        }
        template<typename Generator>
        Point operator()(Generator& gen, typename base::param_type const& p)
        {
            std::vector<std::size_t> skip_list;
            std::vector<domain_point_type> point_cache;
            std::vector<length_type> accumulated_lengths;
            init(skip_list, point_cache, accumulated_lengths, p);
            typedef typename select_most_precise<double, typename coordinate_type<Point>::type>::type sample_type;
            std::uniform_real_distribution<sample_type> dist(0, 1);
            return sample_multi_line<Point>(point_cache, skip_list, accumulated_lengths, dist(gen)*accumulated_lengths.back());
        }
        template<typename Generator>
        Point operator()(Generator& gen)
        {
            typedef typename select_most_precise<double, typename coordinate_type<Point>::type>::type sample_type;
            std::uniform_real_distribution<sample_type> dist(0, 1);
            return sample_multi_line<Point>(point_cache, skip_list, accumulated_lengths, dist(gen)*accumulated_lengths.back());
        }
    private:
        std::vector<std::size_t> skip_list;
        std::vector<domain_point_type> point_cache;
        std::vector<length_type> accumulated_lengths;
    };

    template<typename Point, typename DomainGeometry, typename SingleOrMulti, std::size_t dimension>
    class uniform_point_distribution<Point, DomainGeometry, segment_tag, SingleOrMulti, cartesian_tag, dimension>
        : public uniform_point_distribution_base<Point, DomainGeometry>
    {
        typedef uniform_point_distribution_base<Point, DomainGeometry> base;
        typedef uniform_point_distribution<Point, segment_view<DomainGeometry>,
            linear_tag, SingleOrMulti, cartesian_tag, dimension> delegate_dist;
    public:
        using base::base;
        using base::param;
        uniform_point_distribution(DomainGeometry const& g):base(g),d(segment_view<DomainGeometry>(g)) {}
        void param(typename base::param_type const& p) {
            d = delegate_dist(segment_view<DomainGeometry>(p.domain()));
        }
        template<typename Generator>
        Point operator()(Generator& gen, typename base::param_type const& p)
        {
            delegate_dist d(segment_view<DomainGeometry>(p.domain()));
            return d(gen);
        }
        template<typename Generator>
        Point operator()(Generator& gen)
        {
            return d(gen);
        }
    private:
        delegate_dist d;
    };

}}}} // namespace boost::geometry::random::detail

template<typename Char, typename Traits, typename Point, typename DomainGeometry>
inline std::basic_ostream<Char, Traits>& operator<<
	(
    	std::basic_ostream<Char, Traits> &os,
        boost::geometry::random::detail::uniform_point_distribution_base<Point, DomainGeometry> const& dist
    )
{
	os << boost::geometry::wkt<DomainGeometry>(dist.domain());
    return os;
}

template<typename Char, typename Traits, typename Point, typename DomainGeometry>
inline std::basic_istream<Char, Traits>& operator>>
	(
		std::basic_istream<Char, Traits> &is, 
		boost::geometry::random::detail::uniform_point_distribution_base<Point, DomainGeometry> & dist
	)
{
	std::basic_string<Char, Traits> line;
	std::getline(is, line);
	DomainGeometry g;
	boost::geometry::read_wkt<DomainGeometry>(line, g);
	dist.param( typename boost::geometry::random::detail::uniform_point_distribution_base<Point, DomainGeometry>::param_type(g) );
	return is;
}
#endif // BOOST_GEOMETRY_EXTENSIONS_RANDOM_DETAIL_UNIFORM_POINT_DISTRIBUTION_HPP
