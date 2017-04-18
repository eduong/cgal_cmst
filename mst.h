#ifndef MST_H_
#define MST_H_

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "GraphDefs.h"

struct EdgeWeightComparator {
	bool operator() (SimpleEdge e1, SimpleEdge e2) {
		return e1.weight < e2.weight;
	}
} EWC;

BoostGraph* computeCustomMst(BoostGraph* g, boost::unordered_set<SimpleEdge>* contraintEdgeSet = NULL) {
	// Inherit vertices from g
	BoostGraph* m = new BoostGraph();

	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		Vertex mV = add_vertex(*m);
		(*m)[mV].pt = (*g)[v].pt;
	}

	// Sort edges by weight in heap
	std::vector<SimpleEdge> edgeVec;
	edgeVec.reserve((*g).m_edges.size());

	// Create 1-to-1 SimpleEdge for each Edge from g
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		CGALPoint cgal_u = (*g)[src].pt;
		CGALPoint cgal_v = (*g)[tar].pt;

		SimpleEdge se(cgal_u, cgal_v, src, tar, 0);
		if (contraintEdgeSet == NULL || contraintEdgeSet->count(se) <= 0) {
			se.weight = CGAL::squared_distance(cgal_u, cgal_v);
		}

		edgeVec.push_back(se);
	}

	std::sort(edgeVec.begin(), edgeVec.end(), EWC);

	// Check for connectivity
	std::vector<int> rank(edgeVec.size());
	std::vector<int> parent(edgeVec.size());
	boost::disjoint_sets<int*, int*, boost::find_with_full_path_compression> ds(&rank[0], &parent[0]);
	for (int i = 0; i < rank.size(); i++) {
		rank[i] = i;
		parent[i] = i;
	}

	for (std::vector<SimpleEdge>::iterator it = edgeVec.begin(); it != edgeVec.end(); ++it) {
		SimpleEdge e = *it;
		Vertex u = e.u_idx;
		Vertex v = e.v_idx;
		//std::cout << (*g)[e].weight << " (" << (*g)[u].pt << ") (" << (*g)[v].pt << ")" << std::endl;
		if (ds.find_set(u) != ds.find_set(v)) {
			ds.link(u, v);
			std::pair<Edge, bool> result = add_edge(u, v, *m);
			assert(result.second);
			Edge mE = result.first;
			(*m)[mE].weight = e.weight;
		}
	}

	return m;
}

/**
* Very slow due to EdgeList accessing Edge properties in linear time
*/
BoostGraph* computeBoostMst(BoostGraph* g) {
	std::list<Edge>* mst = new std::list<Edge>();

	boost::kruskal_minimum_spanning_tree(
		*g,
		std::back_inserter(*mst),
		boost::weight_map(get(&EdgeProperties::weight, *g))
		);

	BoostGraph* m = new BoostGraph();

	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		Vertex mV = add_vertex(*m);
		(*m)[mV].pt = (*g)[v].pt;
	}

	for (std::list<Edge>::iterator it = mst->begin(); it != mst->end(); ++it) {
		Edge e = *it;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		std::pair<Edge, bool> result = add_edge(src, tar, *m);
		assert(result.second);
		Edge mE = result.first;
		(*m)[mE].weight = (*g)[e].weight;
	}

	delete mst;
	return m;
}

#endif  // MST_H_