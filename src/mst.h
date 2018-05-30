#ifndef MST_H_
#define MST_H_

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "GraphDefs.h"

struct EdgeWeightComparator2 {
	boost::hash<SimpleEdge> hasher;

	bool operator() (SimpleEdge* e1, SimpleEdge* e2) {
		// The sort order must be deterministic even if weights are the same.
		// The reason for this is: suppose you have a graph with 3 edges, each with equal distance. 
		// A construction of MST is possible with any 2 edges. If sort order is not deterministic
		// picking any 2 of these edges is valid. Ultimately, that changes the outcome of the MST, which leads
		// to non-deterministic results later in the algorithm.
		// Hence, if weights are the same, we sort the edge by hash value, an arbitrary
		// comparison but sufficient for keeping the order deterministic while guarding against u, v index
		// being swapped , i.e. (u, v) is the same edge as (v, u), during the dynamic tree cycle check.
		if (e1->weight == e2->weight) {
			return hasher(*e1) < hasher(*e2);
		}
		return e1->weight < e2->weight;
	}
} EWC2;

std::vector<SimpleEdge*>* sortByWeight(VertexVector* vertices, EdgeVector* edges, boost::unordered_set<SimpleEdge>* contraintEdgeSet = NULL) {
	std::vector<SimpleEdge*>* edgeVec = new std::vector<SimpleEdge*>();
	edgeVec->reserve(edges->size());

	// TODO Simplify this. No longer needed to copy edges.
	// Create 1-to-1 SimpleEdge for each Edge from g. A SimpleEdge has a hash function that maps
	// the 2 endpoints to the same hash value. e.g. Hash value for (1, 2)(9, 10) is the same as (9, 10)(1, 2)
	// Simple Edges also have a weight property that can be sorted on via EdgeWeightComparator
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		CGALPoint* cgal_u = (*vertices)[edge->u];
		CGALPoint* cgal_v = (*vertices)[edge->v];

		SimpleEdge* se = new SimpleEdge(edge->u, edge->v, 0);
		// Recalculate the edge weight for edges not in contraintEdgeSet
		// Edges found in the constraintEdgeSet assume 0 edge weight
		if (contraintEdgeSet == NULL || contraintEdgeSet->count(*se) <= 0) {
			se->weight = CGAL::to_double(CGAL::squared_distance(*cgal_u, *cgal_v));
		}

		edgeVec->push_back(se);
	}

	std::sort(edgeVec->begin(), edgeVec->end(), EWC2);

	return edgeVec;
}

EdgeVector* computeCustomMst(VertexVector* vertices, EdgeVector* edges, boost::unordered_set<SimpleEdge>* contraintEdgeSet = NULL) {
	// Sort edges by weight in heap
	std::vector<SimpleEdge*>* edgeVec = sortByWeight(vertices, edges, contraintEdgeSet);

	// Assign sorted order
	for (int i = 0; i < edgeVec->size(); i++) {
		SimpleEdge* e = (*edgeVec)[i];
		double w = to_double(e->weight);
		e->sortedOrder = i + 1;
	}

	// Init connectivity check data structure
	std::vector<int> rank(edgeVec->size());
	std::vector<int> parent(edgeVec->size());
	boost::disjoint_sets<int*, int*, boost::find_with_full_path_compression> ds(&rank[0], &parent[0]);
	for (int i = 0; i < rank.size(); i++) {
		rank[i] = i;
		parent[i] = i;
	}

	EdgeVector* mst = new EdgeVector();
	mst->reserve(vertices->size() - 1); // There are vertices->size() - 1 edges in an mst

	for (std::vector<SimpleEdge*>::iterator it = edgeVec->begin(); it != edgeVec->end(); ++it) {
		SimpleEdge* e = *it;
		VertexIndex u = e->u;
		VertexIndex v = e->v;

		//std::cout << (*g)[e].weight << " (" << (*g)[u].pt << ") (" << (*g)[v].pt << ")" << std::endl;
		if (ds.find_set(u) != ds.find_set(v)) {
			ds.link(u, v);
			mst->push_back(new SimpleEdge(u, v, 0, e->sortedOrder));
		}

		// Mst complete
		if (mst->size() == vertices->size() - 1) {
			break;
		}
	}

	deleteEdgeVector(edgeVec);

	return mst;
}

/**
* Very slow due to EdgeList accessing Edge properties in linear time
* and the Priority Queue seems to have inefficiencies
*/
//BoostGraph* computeBoostMst(BoostGraph* g) {
//	std::list<Edge>* mst = new std::list<Edge>();
//
//	boost::kruskal_minimum_spanning_tree(
//		*g,
//		std::back_inserter(*mst),
//		boost::weight_map(get(&EdgeProperties::weight, *g))
//		);
//
//	BoostGraph* m = new BoostGraph();
//
//	std::pair<VertexIter, VertexIter> vp;
//	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
//		Vertex v = *vp.first;
//		Vertex mV = add_vertex(*m);
//		(*m)[mV].pt = (*g)[v].pt;
//	}
//
//	for (std::list<Edge>::iterator it = mst->begin(); it != mst->end(); ++it) {
//		Edge e = *it;
//		Vertex src = source(e, *g);
//		Vertex tar = target(e, *g);
//		std::pair<Edge, bool> result = add_edge(src, tar, *m);
//		assert(result.second);
//		Edge mE = result.first;
//		(*m)[mE].weight = (*g)[e].weight;
//	}
//
//	delete mst;
//	return m;
//}

#endif  // MST_H_