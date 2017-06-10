#ifndef GRAPH_DEFS_H_
#define GRAPH_DEFS_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/linear_congruential.hpp>

// CGAL typedefs (2 space)
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 CGALPoint;
typedef K::Segment_2 CGALSegment;
typedef K::Intersect_2 CGALIntersect;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, CGAL::No_intersection_tag> CDT;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face Face_handle;

typedef CGAL::Creator_uniform_2<double, CGALPoint> Creator;

typedef double EdgeWeight;

struct VertexProperties
{
	CGALPoint pt;
};

struct EdgeProperties
{
	EdgeWeight weight;
};

bool operator==(CGALPoint const& p1, CGALPoint const& p2)
{
	return p1.x() == p2.x() && p1.y() == p2.y();
}

namespace boost {
	template <> struct hash < CGALPoint > {
		size_t operator()(CGALPoint const& p) const {
			std::size_t seed = 31;
			boost::hash_combine(seed, p.x());
			boost::hash_combine(seed, p.y());
			return seed;
		}
	};
}

struct SimpleEdge {
	CGALPoint u;
	CGALPoint v;
	size_t u_idx;
	size_t v_idx;
	EdgeWeight weight;

	SimpleEdge(CGALPoint aU, CGALPoint aV, size_t aU_idx, size_t aV_idx, int aWeight) {
		u = aU;
		v = aV;
		u_idx = aU_idx;
		v_idx = aV_idx;
		weight = aWeight;
	}
};

bool operator==(SimpleEdge const& e1, SimpleEdge const& e2)
{
	return (e1.u == e2.u && e1.v == e2.v) || (e1.u == e2.v && e1.v == e2.u);
}

namespace boost {
	template <> struct hash < SimpleEdge > {
		size_t operator()(SimpleEdge const& e) const {
			std::size_t seed = 31;
			boost::hash<CGALPoint> point_hasher;

			// Add hash so edge endpoints {(x1, y1) (x2, y2)}
			// have the same hash as {(x2, y2) (x1, y1)}
			return point_hasher(e.u) + point_hasher(e.v);
		}
	};
}

// Custom data types
typedef std::vector<CGALPoint> VertexVector;
typedef size_t VertexIndex;
typedef std::vector<std::pair<VertexIndex, VertexIndex>*> EdgeVector;
typedef boost::unordered_map<int, Vertex_handle> VertexHandleMap;

struct SimpleEdge2 {
	VertexIndex u_idx;
	VertexIndex v_idx;
	EdgeWeight weight;

	SimpleEdge2(VertexIndex aU_idx, VertexIndex aV_idx, int aWeight) {
		u_idx = aU_idx;
		v_idx = aV_idx;
		weight = aWeight;
	}
};

bool operator==(SimpleEdge2 const& e1, SimpleEdge2 const& e2)
{
	return (e1.u_idx == e2.u_idx && e1.v_idx == e2.v_idx) || (e1.u_idx == e2.v_idx && e1.v_idx == e2.u_idx);
}

namespace boost {
	template <> struct hash < SimpleEdge2 > {
		size_t operator()(SimpleEdge2 const& e) const {
			std::size_t seed = 31;
			boost::hash<VertexIndex> index_hasher;

			// Add hash so edge endpoints {(x1, y1) (x2, y2)}
			// have the same hash as {(x2, y2) (x1, y1)}
			return index_hasher(e.u_idx) + index_hasher(e.v_idx);
		}
	};
}

typedef boost::adjacency_list <
	boost::hash_setS, // OutEdgeList
	boost::vecS, // VertexList
	boost::undirectedS, // Undirected edges
	VertexProperties, // Vertex obj representation
	EdgeProperties, // Edge obj representation
	boost::no_property, // Graph obj representation
	boost::listS > // EdgeList (there is limited documentation on how to change this)
	BoostGraph;

// Boost typdefs
typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<BoostGraph>::edge_descriptor Edge;

typedef boost::graph_traits<BoostGraph>::vertex_iterator VertexIter;
typedef boost::graph_traits<BoostGraph>::edge_iterator EdgeIter;

// Random gen
boost::minstd_rand gen;

#endif  // GRAPH_DEFS_H_