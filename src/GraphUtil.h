#ifndef GRAPH_UTIL_H_
#define GRAPH_UTIL_H_

#include "GraphDefs.h"

boost::unordered_set<SimpleEdge>* createSimpleEdgeSet(BoostGraph* g) {
	boost::unordered_set<SimpleEdge>* contraintEdgeSet = new boost::unordered_set<SimpleEdge>();

	// Map each edge into a hash-able edge
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		CGALPoint cgal_u = (*g)[src].pt;
		CGALPoint cgal_v = (*g)[tar].pt;

		SimpleEdge se(cgal_u, cgal_v, src, tar, 0);
		contraintEdgeSet->insert(se);
	}
	return contraintEdgeSet;
}

boost::unordered_set<SimpleEdge>* createSimpleEdgeSet(VertexVector* vertices, EdgeVector* edges) {
	boost::unordered_set<SimpleEdge>* contraintEdgeSet = new boost::unordered_set<SimpleEdge>();

	// Map each edge into a hash-able edge
	for (int i = 0; i < edges->size(); i++) {
		std::pair<VertexIndex, VertexIndex> edge = (*edges)[i];
		VertexIndex u = edge.first;
		VertexIndex v = edge.second;
		CGALPoint cgal_u = (*vertices)[u];
		CGALPoint cgal_v = (*vertices)[v];

		SimpleEdge se(cgal_u, cgal_v, u, v, 0);
		contraintEdgeSet->insert(se);
	}
	return contraintEdgeSet;
}

BoostGraph* CopyVertices(BoostGraph* src) {
	BoostGraph* tar = new BoostGraph();
	//BoostGraph* tar = new BoostGraph(src->m_vertex_set.size());

	// Inherit all vertices (in same index ordered, i.e. a vertex in g corresponds to the same index and coord as in s)
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*src); vp.first != vp.second; ++vp.first) {
		Vertex v = add_vertex(*tar);
		(*tar)[v].pt = (*src)[*vp.first].pt;
	}
	return tar;
}

#endif  // GRAPH_UTIL_H_