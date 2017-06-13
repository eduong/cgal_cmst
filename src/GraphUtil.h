#ifndef GRAPH_UTIL_H_
#define GRAPH_UTIL_H_

#include "GraphDefs.h"

boost::unordered_set<SimpleEdge>* createSimpleEdgeSet(EdgeVector* edges) {
	boost::unordered_set<SimpleEdge>* contraintEdgeSet = new boost::unordered_set<SimpleEdge>();

	// Map each edge into a hash-able edge
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		contraintEdgeSet->insert(*edge);
	}
	return contraintEdgeSet;
}

BoostGraph* CopyVertices(BoostGraph* src) {
	BoostGraph* tar = new BoostGraph();

	// Inherit all vertices (in same index ordered, i.e. a vertex in g corresponds to the same index and coord as in s)
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*src); vp.first != vp.second; ++vp.first) {
		Vertex v = add_vertex(*tar);
		(*tar)[v].pt = (*src)[*vp.first].pt;
	}
	return tar;
}

#endif  // GRAPH_UTIL_H_