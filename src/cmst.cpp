#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Segment_2.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/linear_congruential.hpp>

#include "Forest.h"

#include <fstream>
#include <string>

// CGAL typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 CGALPoint;
typedef K::Segment_2 CGALSegment;

typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, CGAL::No_intersection_tag> CDT;

struct VertexProperties
{
	CGALPoint pt;
};

struct EdgeProperties
{
	int weight;
};

typedef boost::adjacency_list < boost::vecS, boost::vecS,
	boost::undirectedS,
	VertexProperties,
	EdgeProperties >
	BoostGraph;

// Boost typdefs
typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<BoostGraph>::edge_descriptor Edge;

typedef boost::graph_traits<BoostGraph>::vertex_iterator VertexIter;
typedef boost::graph_traits<BoostGraph>::edge_iterator EdgeIter;

boost::minstd_rand gen;

CDT* computeCdt(std::vector<CGALPoint>* vertices, std::list<std::pair<int, int>>* edges);
BoostGraph* convertCdtToGraph(std::vector<CGALPoint>* vertices, CDT* cdt, bool assignZeroWeightToConstraints);
void printGraph(const char* title, BoostGraph* g);

BoostGraph* CreateRandomPlaneForest(int numVertices, int bounds, double edgeProbability) {
	std::vector<CGALPoint>* vertices = new std::vector<CGALPoint>();

	// Define bounds random gen
	boost::uniform_int<> boundsRange(0, bounds);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> boundsDice(gen, boundsRange);
	boundsDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	// To maintain unique x, y
	std::map<int, std::map<int, int>> xBounds;

	for (int i = 0; i < numVertices; i++) {
		int x = boundsDice();
		int y = boundsDice();
		if (xBounds.count(x)) { // x exists
			while (xBounds[x].count(y)) { // y exists, reroll (TODO: range check to prevent infinite loop)
				y = boundsDice();
			}
			xBounds[x].insert(std::pair<int, int>(y, NULL));
		}
		else {
			std::map<int, int> yBounds;
			yBounds.insert(std::pair<int, int>(y, NULL));
			xBounds.insert(std::pair<int, std::map<int, int>>(x, yBounds));
		}
		vertices->push_back(CGALPoint(x, y));
	}

	// Compute DT (no constraints)
	CDT* cdt = computeCdt(vertices, NULL);

	BoostGraph* g = convertCdtToGraph(vertices, cdt, false);
	//printGraph("Random Forest", g);
	vertices->clear();
	delete vertices;

	// Removal of edges on vecS graphs invalidates the iterator. It is easier to through them in a list
	std::list<Edge> edges;
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		edges.push_back(e);
	}

	// Define edge random gen
	boost::uniform_real<> edgeRange(0, 1);
	boost::variate_generator<boost::minstd_rand, boost::uniform_real<>> edgeDice(gen, edgeRange);
	edgeDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	int count = 0;
	int removed = 0;
	for (std::list<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) {
		Edge e = *it;
		// E.g. edgeDice is 0.6 and edgeProbability = 0.35 then the edge will be removed
		double roll = edgeDice();
		if (roll > edgeProbability) {
			remove_edge(e, (*g));
			removed++;
		}
		count++;
	}

	std::cout << "Edge count before removal " << count << std::endl;
	std::cout << "Edge count after removal " << (count - removed) << std::endl;
	std::cout << "Ratio " << (double(count - removed) / (double)count) << std::endl;

	return g;
}

void printGraph(const char* title, BoostGraph* g) {
	std::cout << std::endl << "=== " << title << std::endl;

	// Iterate through the vertices and print them out
	std::cout << "Vertices: " << std::endl;
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		std::cout << "index: " << v << " (" << (*g)[v].pt << ")" << std::endl;
	}

	// Iterate through the edges and print them out
	std::cout << std::endl << "Edges: " << std::endl;
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		std::cout << (*g)[e].weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
	}
	std::cout << std::endl;
}

void printGraph(const char* title, std::vector<CGALPoint>* vertices, std::list<std::pair<int, int>>* edges) {
	std::cout << std::endl << "=== " << title << std::endl;

	// Insert vertices
	std::cout << "Vertices: " << std::endl;
	for (std::vector<CGALPoint>::iterator it = vertices->begin(); it != vertices->end(); ++it) {
		std::cout << "(" << *it << ")" << std::endl;
	}

	std::cout << std::endl << "Edges: " << std::endl;
	// Insert constraint edges
	for (std::list<std::pair<int, int>>::iterator it = edges->begin(); it != edges->end(); ++it) {
		std::cout << "(" << (*vertices)[it->first] << ") (" << (*vertices)[it->second] << ")" << std::endl;
	}
}

CDT* computeCdt(std::vector<CGALPoint>* vertices, std::list<std::pair<int, int>>* edges) {
	CDT* cdt = new CDT();

	// Insert vertices
	for (std::vector<CGALPoint>::iterator it = vertices->begin(); it != vertices->end(); ++it) {
		cdt->insert(*it);
	}

	// Insert constraint edges
	if (edges) {
		for (std::list<std::pair<int, int>>::iterator it = edges->begin(); it != edges->end(); ++it) {
			cdt->insert_constraint((*vertices)[it->first], (*vertices)[it->second]);
		}
	}

	assert(cdt->is_valid());
	return cdt;
}

CDT* computeCdt(BoostGraph* g) {
	CDT* cdt = new CDT();

	// Insert vertices
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		cdt->insert((*g)[v].pt);
	}

	// Insert constraint edges
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		cdt->insert_constraint((*g)[src].pt, (*g)[tar].pt);
	}

	assert(cdt->is_valid());
	return cdt;
}

void printCdtInfo(CDT* cdt) {
	std::cout << std::endl << "=== Contraint Delaunay Triangulation" << std::endl;
	int constraintCount = 0, unconstraintCount = 0;
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge e = *eit;
		if (cdt->is_constrained(e)) {
			constraintCount++;
			std::cout << "C ";
		}
		else {
			unconstraintCount++;
			std::cout << "U ";
		}

		CGALSegment s = cdt->segment(e);
		const CGALPoint& p1 = s.point(0);
		const CGALPoint& p2 = s.point(1);
		std::cout << "(" << p1 << ") (" << p2 << ")" << std::endl;
	}
	std::cout << "# constrained edges " << constraintCount << std::endl;
	std::cout << "# unconstrained edges " << unconstraintCount << std::endl << std::endl;

	for (CDT::Vertex_iterator vit = cdt->vertices_begin(); vit != cdt->vertices_end(); ++vit) {
		CDT::Vertex v = *vit;
		std::cout << "(" << v << ")" << std::endl;
	}
}

BoostGraph* convertCdtToGraph(BoostGraph* g, CDT* cdt, bool assignZeroWeightToConstraints) {
	BoostGraph* bg_cdt = new BoostGraph();

	// Add vertices to graph
	std::map<CGALPoint, Vertex> verticesMap;
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		Vertex bg_v = add_vertex(*bg_cdt);
		(*bg_cdt)[bg_v].pt = (*g)[v].pt;
		verticesMap.insert(std::pair<CGALPoint, Vertex>((*g)[v].pt, bg_v));
	}

	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		Vertex u = verticesMap[cgal_u];
		Vertex v = verticesMap[cgal_v];
		std::pair<Edge, bool> result = add_edge(u, v, *bg_cdt);
		assert(result.second);
		Edge e = result.first;
		if (assignZeroWeightToConstraints && cdt->is_constrained(cgal_e)) {
			(*bg_cdt)[e].weight = 0;
		}
		else {
			(*bg_cdt)[e].weight = CGAL::squared_distance(cgal_u, cgal_v);
		}
	}

	return bg_cdt;
}

BoostGraph* convertCdtToGraph(std::vector<CGALPoint>* vertices, CDT* cdt, bool assignZeroWeightToConstraints) {
	BoostGraph* g = new BoostGraph();

	// Add vertices to graph
	std::map<CGALPoint, Vertex> verticesMap;
	for (std::vector<CGALPoint>::iterator it = vertices->begin(); it != vertices->end(); ++it) {
		Vertex v = add_vertex(*g);
		(*g)[v].pt = *it;
		verticesMap.insert(std::pair<CGALPoint, Vertex>(*it, v));
	}

	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		Vertex u = verticesMap[cgal_u];
		Vertex v = verticesMap[cgal_v];
		std::pair<Edge, bool> result = add_edge(u, v, *g);
		assert(result.second);
		Edge e = result.first;
		if (assignZeroWeightToConstraints && cdt->is_constrained(cgal_e)) {
			(*g)[e].weight = 0;
		}
		else {
			(*g)[e].weight = CGAL::squared_distance(cgal_u, cgal_v);
		}
	}

	return g;
}

BoostGraph* computeMst(BoostGraph* g) {
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

	return m;
}

struct tree_visitor : boost::default_bfs_visitor {
	Forest* forest;
	tree_visitor(Forest* f) {
		forest = f;
	}

	void tree_edge(const Edge &e, const BoostGraph &g) const {
		Vertex src = source(e, g);
		Vertex tar = target(e, g);
		int weight = CGAL::squared_distance(g[src].pt, g[tar].pt); // Re-calculate weights since constraint edges have 0 weight
		forest->Link(tar, src);
		forest->SetCost(tar, weight); // First param to Link() is always the leafmost node
		std::cout << "Tree edge: w:" << weight << " : " << tar << " (" << g[tar].pt << ") " << src << " (" << g[src].pt << ")" << std::endl;
	}
};

Forest* createLinkCutTree(BoostGraph* g) {
	Forest* f = new Forest();
	f->Initialize(g->m_vertices.size());

	tree_visitor vis(f);

	// Arbitrary root vertex
	Vertex s = *(boost::vertices(*g).first);
	boost::breadth_first_search(*g, s, boost::visitor(vis));

	return f;
}

BoostGraph* checkCycles(Forest* f, BoostGraph* g) {
	BoostGraph* s = new BoostGraph();
	//f->ZeroUnitializedCosts();

	// Inherit all vertices (in same index ordered, i.e. a vertex in g corresponds to the same index and coord as in s)
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		Vertex sV = add_vertex(*s);
		(*s)[sV].pt = (*g)[v].pt;
	}

	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		int we = (*g)[e].weight;
		Vertex u = source(e, *g);
		Vertex v = target(e, *g);

		if (!f->SameEdge(u, v) && f->SameTree(u, v)) {
			Forest::NodeId lca = f->LCA(u, v);
			int lcaWeight = f->GetCost(lca); // LCA weight needs to be ignored (zero weight) because all weights are stored in the leafmost node, i.e. the weight on node lca represents edge (lca, parent lca)
			f->SetCost(lca, 0);
			Forest::NodeId parentOfLca = f->Cut(lca);

			Forest::NodeId max_u = f->FindMax(u);
			Forest::Cost cost_u = f->GetCost(max_u);
			while (we < cost_u) {
				f->SetCost(max_u, 0);

				std::pair<Edge, bool> result = add_edge(max_u, f->FindParent(max_u), *s);
				assert(result.second);
				Edge e = result.first;
				Vertex src = source(e, *s);
				Vertex tar = target(e, *s);
				//std::cout << "(" << (*s)[src].pt << ") (" << (*s)[tar].pt << ")" << std::endl;

				max_u = f->FindMax(u);
				cost_u = f->GetCost(max_u);
			}

			Forest::NodeId max_v = f->FindMax(v);
			Forest::Cost cost_v = f->GetCost(max_v);
			while (we < cost_v) {
				f->SetCost(max_v, 0);

				std::pair<Edge, bool> result = add_edge(max_v, f->FindParent(max_v), *s);
				assert(result.second);
				Edge e = result.first;
				Vertex src = source(e, *s);
				Vertex tar = target(e, *s);
				//std::cout << "(" << (*s)[src].pt << ") (" << (*s)[tar].pt << ")" << std::endl;

				max_v = f->FindMax(v);
				cost_v = f->GetCost(max_v);
			}

			f->Link(lca, parentOfLca);
			f->SetCost(lca, lcaWeight);
		}
	}
	return s;
}

BoostGraph* computeCmst(BoostGraph* f) {
	// Input: plane forest F = (V, E)
	// Output: minimum set S ⊆ E of constraints such that F ⊆ CMST(V, S)
	printGraph("Input graph F = (V, E)", f);

	// Compute CDT(F)
	CDT* cdt = computeCdt(f);
	printCdtInfo(cdt);

	// Create graph representations for CDT(F)
	BoostGraph* g1 = convertCdtToGraph(f, cdt, false);
	printGraph("Contraint Delaunay Triangulation (as BoostGraph)", g1);

	// Create graph representations for CDT◦(F)
	BoostGraph* g2 = convertCdtToGraph(f, cdt, true);
	printGraph("CDT◦(F)", g2);

	// Compute T' = MST(CDT(F))
	BoostGraph* mst = computeMst(g1);
	printGraph("T' = MST(CDT(F))", mst);

	// Compute CMST(F) = MST(CDT◦(F))
	BoostGraph* mst2 = computeMst(g2);
	printGraph("CMST(F) = MST(CDT◦(F))", mst2);

	// Create DynamicTree for CMST(F)
	Forest* forest = createLinkCutTree(mst2);

	// Check for cycles and construct constraint set s
	BoostGraph* s = checkCycles(forest, mst);
	printGraph("Contraint set S", s);

	delete forest;
	
	mst2->clear();
	delete mst2;

	mst->clear();
	delete mst;

	g2->clear();
	delete g2;
	
	g1->clear();
	delete g1;

	cdt->clear();
	delete cdt;

	return s;
}

int main(int argc, char* argv[]) {
	BoostGraph* F = CreateRandomPlaneForest(10, 10, 0.20);
	printGraph("Input plane forest F = (V, E)", F);
	
	BoostGraph* S = computeCmst(F);

	// Validate F ⊆ CMST(V, S)
	CDT* cdtF = computeCdt(F);
	CDT* cdtS = computeCdt(S);
	BoostGraph* bg_cdtF = convertCdtToGraph(F, cdtF, true);
	BoostGraph* bg_cdtS = convertCdtToGraph(S, cdtS, true);
	BoostGraph* mstF = computeMst(bg_cdtF);
	BoostGraph* mstS = computeMst(bg_cdtS);
	printGraph("MST(F)", mstF);
	printGraph("MST(S)", mstS);

	delete mstS;
	delete mstF;
	delete bg_cdtS;
	delete bg_cdtF;
	delete cdtS;
	delete cdtF;
	delete S;

	delete F;

	return 0;
}