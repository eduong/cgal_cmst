﻿#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/intersections.h>
#include <CGAL/intersections_d.h>

#include <boost/chrono.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/heap/priority_queue.hpp>

#include "Forest.h"
#include "GraphDefs.h"
#include "GraphUtil.h"
#include "mst.h"
#include "timer.h"
#include "verify.h"

#include <fstream>
#include <string>

#define SHOW_DEBUG true

/**
* Naive linear time intersection
* Returns true if edge (u, v) intersects an edge in g, otherwise false
**/
bool DoesIntersect(BoostGraph* g, Vertex u, Vertex v) {
	CGALPoint uPt = (*g)[u].pt;
	CGALPoint vPt = (*g)[v].pt;
	CGALSegment segUV(uPt, vPt);

	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		CGALSegment seg((*g)[src].pt, (*g)[tar].pt);

		CGAL::cpp11::result_of<CGALIntersect(CGALSegment, CGALSegment)>::type result = intersection(seg, segUV);
		if (result) {
			if (const CGALSegment* s = boost::get<CGALSegment>(&*result)) {
				//std::cout << *s << std::endl;
				return true;
			}
			else if (const CGALPoint* p = boost::get<CGALPoint >(&*result)) {
				//std::cout << " i " << *p;
				// Ignore intersection at segment endpoints
				if (*p != uPt && *p != vPt) {
					return true;
				}
			}
		}
	}
	return false;
}

BoostGraph* CreateRandomPlaneForest(int numVertices, int bounds, int edgeRolls, double edgeProbability) {
	BoostGraph* g = new BoostGraph();

	// Define bounds random gen
	boost::uniform_int<> boundsRange(0, bounds);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> boundsDice(gen, boundsRange);
	boundsDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	// To maintain unique x, y
	boost::unordered_map<int, boost::unordered_map<int, int>> xBounds;

	// Generate vertices with random coordinated within bounds
	for (int i = 0; i < numVertices; i++) {
		int x = boundsDice();
		int y = boundsDice();
		if (xBounds.count(x)) { // x exists
			while (xBounds[x].count(y)) { // y exists, reroll (TODO: range check to prevent infinite loop)
				y = boundsDice();
			}
			xBounds[x].emplace(y, NULL);
		}
		else {
			boost::unordered_map<int, int> yBounds;
			yBounds.emplace(y, NULL);
			xBounds.emplace(x, yBounds);
		}
		Vertex v = add_vertex(*g);
		(*g)[v].pt = CGALPoint(x, y);
	}

	// Define edge random gen
	boost::uniform_real<> edgeRange(0, 1);
	boost::variate_generator<boost::minstd_rand, boost::uniform_real<>> edgeDice(gen, edgeRange);
	edgeDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	boost::uniform_int<> vertexRange(0, numVertices - 1);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> vertexDice(gen, vertexRange);
	vertexDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	std::vector<int> rank(numVertices);
	std::vector<int> parent(numVertices);
	boost::disjoint_sets<int*, int*, boost::find_with_full_path_compression> ds(&rank[0], &parent[0]);
	for (int i = 0; i < rank.size(); i++) {
		rank[i] = i;
		parent[i] = i;
	}

	// Select random vertices u, v for edgeRolls number of times
	// An edge connects u, v:
	//		1. u != v
	//		2. roll <= edgeProbability
	//		3. adding edge(u, v) does not create a cycle
	//		4. edge(u, v) does not intersect any other edge
	for (int i = 0; i < edgeRolls; i++) {
		Vertex u = vertexDice();
		Vertex v = vertexDice();
		double d = edgeDice();
		if (u != v
			&& d <= edgeProbability
			&& ds.find_set(u) != ds.find_set(v)
			&& !DoesIntersect(g, u, v)) {

			// Add edge(u, v)
			std::pair<Edge, bool> result = add_edge(u, v, *g);
			assert(result.second);
			ds.link(u, v);
			//std::cout << " - added";
		}
	}

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
		//std::cout << edgeWeightMap->at(e).weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
		std::cout << (*g)[e].weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
	}
	std::cout << std::endl;
}

CDT* computeCdt(BoostGraph* g) {
	CDT* cdt = new CDT();

	// Insert vertices
	boost::unordered_map<CGALPoint, Vertex_handle> vertexHandles;
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		CGALPoint pt = (*g)[v].pt;
		Vertex_handle vHandle = cdt->insert(pt);
		vertexHandles.emplace(pt, vHandle);
	}

	// Insert constraint edges
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex u = source(e, *g);
		Vertex v = target(e, *g);
		Vertex_handle uH = vertexHandles[(*g)[u].pt];
		Vertex_handle vH = vertexHandles[(*g)[v].pt];
		cdt->insert_constraint(uH, vH);
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

	/*for (CDT::Vertex_iterator vit = cdt->vertices_begin(); vit != cdt->vertices_end(); ++vit) {
		CDT::Vertex v = *vit;
		std::cout << "(" << v << ")" << std::endl;
		}*/
}

BoostGraph* convertCdtToGraph(BoostGraph* g, CDT* cdt) {
	BoostGraph* bg_cdt = new BoostGraph();

	// Add vertices to graph
	boost::unordered_map<CGALPoint, Vertex> verticesMap;
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		Vertex bg_v = add_vertex(*bg_cdt);
		(*bg_cdt)[bg_v].pt = (*g)[v].pt;
		verticesMap.emplace((*g)[v].pt, bg_v);
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		Vertex u = verticesMap[cgal_u];
		Vertex v = verticesMap[cgal_v];
		std::pair<Edge, bool> result = add_edge(u, v, *bg_cdt);
		assert(result.second);
	}

	return bg_cdt;
}

struct tree_visitor : boost::default_bfs_visitor {
public:
	Forest* forest;

	tree_visitor(Forest* f) : forest(f) {}

	void tree_edge(const Edge &e, const BoostGraph &g) const {
		Vertex src = source(e, g);
		Vertex tar = target(e, g);
		EdgeWeight weight = CGAL::squared_distance(g[src].pt, g[tar].pt); // Re-calculate weights since constraint edges have 0 weight
		forest->Link(tar, src);
		forest->SetCost(tar, weight); // First param to Link() is always the leafmost node
		//std::cout << "Tree edge: w:" << weight << " : " << tar << " (" << g[tar].pt << ") " << src << " (" << g[src].pt << ")" << std::endl;
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
	BoostGraph* s = InheritVertices(g);

	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		EdgeWeight we = (*g)[e].weight;
		Vertex u = source(e, *g);
		Vertex v = target(e, *g);

		if (!f->SameEdge(u, v) && f->SameTree(u, v)) {
			Forest::NodeId lca = f->LCA(u, v);
			EdgeWeight lcaWeight = f->GetCost(lca); // LCA weight needs to be ignored (zero weight) because all weights are stored in the leafmost node, i.e. weight on node lca belongs to edge (lca, parent lca)
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

void computeCmst(BoostGraph* F, BoostGraph** TPrime, BoostGraph** Cmst, BoostGraph** S) {
	// Input: plane forest F = (V, E)
	// Output: minimum set S ⊆ E of constraints such that F ⊆ CMST(V, S)
	if (SHOW_DEBUG) { printGraph("Input graph F = (V, E)", F); }

	boost::unordered_set<SimpleEdge>* contraintEdgeSet = createSimpleEdgeSet(F);

	boost::chrono::high_resolution_clock::time_point start;
	boost::chrono::high_resolution_clock::time_point end;
	boost::chrono::milliseconds duration(0);
	boost::chrono::milliseconds total(0);

	// Compute CDT(F)
	start = boost::chrono::high_resolution_clock::now();
	CDT* cdt = computeCdt(F);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("CDT(F)", duration);
	if (SHOW_DEBUG) { printCdtInfo(cdt); }

	// Create boost graph representations for CDT(F)
	start = boost::chrono::high_resolution_clock::now();
	BoostGraph* g1 = convertCdtToGraph(F, cdt);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Graph conversion CDT(F)", duration);
	if (SHOW_DEBUG) { printGraph("Contraint Delaunay Triangulation (as BoostGraph)", g1); }

	// Compute T' = MST(CDT(F)) (Best case MST)
	start = boost::chrono::high_resolution_clock::now();
	*TPrime = computeCustomMst(g1);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("T' = MST(CDT(F))", duration);
	if (SHOW_DEBUG) { printGraph("T' = MST(CDT(F))", *TPrime); }

	// Compute CMST(F) = MST(CDT◦(F))
	start = boost::chrono::high_resolution_clock::now();
	//*Cmst = computeCustomMst(g2);
	*Cmst = computeCustomMst(g1, contraintEdgeSet);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("CMST(F) = MST(CDT◦(F))", duration);
	if (SHOW_DEBUG) { printGraph("CMST(F) = MST(CDT◦(F))", *Cmst); }

	// Create DynamicTree for CMST(F)
	start = boost::chrono::high_resolution_clock::now();
	Forest* forest = createLinkCutTree(*Cmst);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Dynamic Tree construction", duration);

	// Check for cycles and construct constraint set s
	start = boost::chrono::high_resolution_clock::now();
	*S = checkCycles(forest, *TPrime);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Check cycles", duration);
	if (SHOW_DEBUG) { printGraph("Contraint set S", *S); }

	printDuration("Total", total);

	delete forest;
	delete contraintEdgeSet;
	delete g1;
	delete cdt;
}

bool containsEdge(BoostGraph* g, boost::unordered_set<SimpleEdge>* edgeSet, Edge e) {
	Vertex u = source(e, *g);
	Vertex v = target(e, *g);

	CGALPoint cgal_u = (*g)[u].pt;
	CGALPoint cgal_v = (*g)[v].pt;

	SimpleEdge se(cgal_u, cgal_v, u, v, 0);
	return edgeSet->count(se) > 0;
}

// True if A a subgraph of B
bool isSubgraph(BoostGraph* a, BoostGraph* b) {
	boost::unordered_set<SimpleEdge>* bEdgeSet = createSimpleEdgeSet(b);

	// Iterate through the edges
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*a); ei != ei_end; ++ei) {
		Edge e = *ei;
		if (!containsEdge(b, bEdgeSet, e)) {
			return false;
		}
	}

	return true;
}

BoostGraph* graphOmitEdge(BoostGraph* S, int omitIndex) {
	BoostGraph* newS = new BoostGraph();

	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*S); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		Vertex nV = add_vertex(*newS);
		(*newS)[nV].pt = (*S)[v].pt;
	}

	int index = 0;
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*S); ei != ei_end; ++ei) {
		if (index++ == omitIndex) {
			continue;
		}

		Edge e = *ei;
		Vertex src = source(e, *S);
		Vertex tar = target(e, *S);

		std::pair<Edge, bool> result = add_edge(src, tar, *newS);
		assert(result.second);
		Edge nE = result.first;
		(*newS)[nE].weight = (*S)[e].weight;
	}

	return newS;
}

bool isCmstSubgraph(BoostGraph* F, BoostGraph* S) {
	CDT* cdtS = computeCdt(S);
	BoostGraph* bg_cdtS = convertCdtToGraph(S, cdtS);
	BoostGraph* cmstS = computeCustomMst(bg_cdtS, createSimpleEdgeSet(S));
	if (SHOW_DEBUG) {
		printGraph("CMST(V, S)", cmstS);
	}

	bool res = isSubgraph(F, cmstS);

	delete cmstS;
	delete bg_cdtS;
	delete cdtS;

	return res;
}

bool isMinimal(BoostGraph* F, BoostGraph* S) {
	// Validate 
	// S ⊆ E s.t. F ⊆ CMST(V, S). Note: CMST(G) = MST(CVG(G)) = MST(CDT◦(G)) (where CDT◦(G) = CDT(G) when all edges in G have 0 weight)
	// Notice: CMST(V, S) is a spanning graph, that is there is at least 1 edge that connects every vertex. So F, the constraint set should be a subset of CMST(V, S)
	// In other words, we want to find the smallest subset S of edges of F such that
	// CMST(F) is equal to CMST(V, S), although the weights of the two trees may
	// be different.
	if (!isCmstSubgraph(F, S)) {
		return false;
	}

	// Removal of any edge of S should result in F !⊆ CMST(V, S)
	for (int i = 0; i < S->m_edges.size(); i++) {
		BoostGraph* omittedS = graphOmitEdge(S, i);
		bool sub = isCmstSubgraph(F, omittedS);
		delete omittedS;

		if (sub) {
			return false;
		}
	}
	return true;
}

bool isContainedIn(BoostGraph* F, BoostGraph* S, BoostGraph* TPrime) {
	boost::unordered_set<SimpleEdge>* sEdgeSet = createSimpleEdgeSet(S);
	boost::unordered_set<SimpleEdge>* tPrimeEdgeSet = createSimpleEdgeSet(TPrime);

	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*F); ei != ei_end; ++ei) {
		Edge e = *ei;
		if (containsEdge(S, sEdgeSet, e))
			continue;
		else if (containsEdge(TPrime, tPrimeEdgeSet, e))
			continue;
		else
			return false;
	}

	return true;
}

// TODO:
// CDT constraint list is sometimes different from F's edge list, is this an CGAL bug?
// Traffic data

int main(int argc, char* argv[]) {

	for (;;) {
		BoostGraph* F = CreateRandomPlaneForest(10, 100000000, 5000, 0.50);
		std::cout << "Random plane forest create with edge count: " << F->m_edges.size() << std::endl;

		BoostGraph* TPrime = NULL;
		BoostGraph* Cmst = NULL;
		BoostGraph* S = NULL;

		computeCmst(F, &TPrime, &Cmst, &S);

		printGraph("Constraint Set S", S);

		// S ⊆ E s.t. F ⊆ CMST(V, S)
		assert(isMinimal(F, S));

		// For each e in F, if e not in S, then e in T'
		assert(isContainedIn(F, S, TPrime));

		delete S;
		delete Cmst;
		delete TPrime;
		delete F;
	}

	return 0;
}