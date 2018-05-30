#include <boost/chrono.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/heap/priority_queue.hpp>

#include <math.h>

#include "Forest.h"
#include "GraphDefs.h"
#include "GraphUtil.h"
#include "GraphParse.h"
#include "GraphGen.h"
#include "Mst.h"
#include "Timer.h"

#include <string>

#define SHOW_DEBUG false
#define SHOW_INPUT_SYNTAX false

void printGraph(const char* title, BoostGraph* g, bool printVertices = false) {
	std::cout << std::endl << "=== " << title << std::endl;

	if (printVertices) {
		// Iterate through the vertices and print them out
		std::cout << "Vertices: " << std::endl;
		std::pair<VertexIter, VertexIter> vp;
		for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
			Vertex v = *vp.first;
			std::cout << "index: " << v << " (" << (*g)[v].pt.x() << ", " << (*g)[v].pt.y() << ")" << std::endl;
		}

		if (SHOW_INPUT_SYNTAX) {
			for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
				Vertex v = *vp.first;
				std::cout << "Vertex v" << v << " = add_vertex(*F); (*F)[v" << v << "].pt = CGALPoint(" << (*g)[v].pt.x() << ", " << (*g)[v].pt.y() << ");" << std::endl;
			}

			std::pair<EdgeIter, EdgeIter> ep;
			EdgeIter ei, ei_end;
			for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
				Edge e = *ei;
				Vertex src = source(e, *g);
				Vertex tar = target(e, *g);
				std::cout << "add_edge(v" << src << ", v" << tar << ", *F);" << std::endl;
			}
		}
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

void printGraph(const char* title, VertexVector* vertices, EdgeVector* edges, bool printVertices = false) {
	std::cout << std::endl << "=== " << title << std::endl;

	if (printVertices) {
		// Iterate through the vertices and print them out
		for (int i = 0; i < vertices->size(); i++) {
			CGALPoint* v = (*vertices)[i];
			std::cout << "index: " << v << " (" << (*v) << ")" << std::endl;
		}
	}

	// Iterate through the edges and print them out
	std::cout << std::endl << "Edges: " << std::endl;
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* e = (*edges)[i];
		CGALPoint* src = (*vertices)[e->u];
		CGALPoint* tar = (*vertices)[e->v];
		//std::cout << edgeWeightMap->at(e).weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
		std::cout << e->weight << " (" << (*src) << ") (" << (*tar) << ")" << std::endl;
	}

	std::cout << std::endl;
}

//#include <iostream>
//#include <fstream>

void printGDFGraph(const char* fileName, VertexVector* vertices, EdgeVector* edges) {
	std::ofstream myfile;
	myfile.open(fileName, std::ios::out | std::ios::in);

	myfile << "nodedef> name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR" << std::endl;
	//std::cout << "nodedef> name VARCHAR,label VARCHAR,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR" << std::endl;

	// Iterate through the vertices and print them out
	boost::unordered_map<CGALPoint*, VertexIndex> vertexHandles;
	for (int i = 0; i < vertices->size(); i++) {
		CGALPoint* v = (*vertices)[i];
		myfile << i << ",,10.0,10.0," << (*v).x() << "," << (*v).y() << ",'153,153,153'" << std::endl;
		//std::cout << i << ",,10.0,10.0," << (*v).x() << "," << (*v).y() << ",'153,153,153'" << std::endl;
		vertexHandles.emplace(v, i);
	}

	myfile << "edgedef> node1,node2,weight DOUBLE,directed BOOLEAN,color VARCHAR" << std::endl;
	//std::cout << "edgedef> node1,node2,weight DOUBLE,directed BOOLEAN,color VARCHAR" << std::endl;

	// Iterate through the edges and print them out
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* e = (*edges)[i];
		CGALPoint* src = (*vertices)[e->u];
		CGALPoint* tar = (*vertices)[e->v];
		VertexIndex srcInd = vertexHandles[src];
		VertexIndex tarInd = vertexHandles[tar];
		//EdgeWeight weight = sqrt(CGAL::squared_distance(*src, *tar));
		EdgeWeight weight = 1.0;
		myfile << srcInd << "," << tarInd << "," << weight << ",false,'128,128,128'" << std::endl;
		//std::cout << srcInd << "," << tarInd << "," << weight << ",false,'128,128,128'" << std::endl;
	}

	myfile.close();
	//std::cout << std::endl;
}

CDT* computeCdt(VertexVector* vertices, EdgeVector* edges) {
	CDT* cdt = new CDT();

	boost::unordered_map<VertexIndex, Vertex_handle> vertexHandles;
	for (int i = 0; i < vertices->size(); i++) {
		CGALPoint* pt = (*vertices)[i];
		Vertex_handle vHandle = cdt->insert(*pt);
		vertexHandles.emplace(i, vHandle);
	}

	// Insert constraint edges
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		VertexIndex u = edge->u;
		VertexIndex v = edge->v;
		Vertex_handle uH = vertexHandles[u];
		Vertex_handle vH = vertexHandles[v];
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

// Given 3 colinear points, if the 2 furthest points are constrained, how is this delaunay triangulated?
// It seems the outcome is implementation dependent. CGAL will split the constraint into 2 edges.
// Other implementations might allow for overlapping, collinear edges. Either way, the current plan
// is to assume that an input forest, F, with collinear edges will be replaced by smaller constraint edges.
EdgeVector* newConstraintSetFromCdt(CDT* cdt, VertexVector* originalVertices) {
	EdgeVector* newEdgeVector = new EdgeVector();

	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < originalVertices->size(); i++) {
		vertexIndex[*(*originalVertices)[i]] = i;
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;

		if (cdt->is_constrained(cgal_e)) {
			// Assumes point coord are unique
			CGALSegment segement = cdt->segment(cgal_e);
			CGALPoint cgal_u = segement.point(0);
			CGALPoint cgal_v = segement.point(1);
			VertexIndex u = vertexIndex[cgal_u];
			VertexIndex v = vertexIndex[cgal_v];

			SimpleEdge* edge = new SimpleEdge(u, v, 0);
			newEdgeVector->push_back(edge);
		}
	}
	return newEdgeVector;
}

EdgeVector* convertCdtToGraph(VertexVector* vertices, CDT* cdt) {
	EdgeVector* edgeVec = new EdgeVector();

	// Map CGALPoint -> VertexIndex
	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < vertices->size(); i++) {
		vertexIndex[*(*vertices)[i]] = i;
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		VertexIndex u = vertexIndex[cgal_u];
		VertexIndex v = vertexIndex[cgal_v];

		SimpleEdge* edge = new SimpleEdge(u, v, 0);
		edgeVec->push_back(edge);
	}

	return edgeVec;
}

BoostGraph* convertVectorToGraph(VertexVector* vertices, EdgeVector* edges, boost::unordered_map<Edge, SimpleEdge*>** boostToSimpleEdgeMap, boost::unordered_map<SimpleEdge, SimpleEdge*>* relativeWeightSet) {
	BoostGraph* g = new BoostGraph();

	// Add vertices to graph
	boost::unordered_map<VertexIndex, Vertex> verticesMap;
	for (int i = 0; i < vertices->size(); i++) {
		CGALPoint* pt = (*vertices)[i];
		Vertex gv = add_vertex(*g);
		(*g)[gv].pt = *pt;
		verticesMap.emplace(i, gv);
	}

	// Map BoostEdge -> SimpleEdge
	(*boostToSimpleEdgeMap) = new boost::unordered_map<Edge, SimpleEdge*>();

	// Add edges to graph
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		Vertex u = verticesMap[edge->u];
		Vertex v = verticesMap[edge->v];
		std::pair<Edge, bool> result = add_edge(u, v, *g);
		assert(result.second);

		SimpleEdge* cdtEdge = (*relativeWeightSet)[(*edge)];
		(**boostToSimpleEdgeMap)[result.first] = cdtEdge;
	}

	return g;
}

struct tree_visitor : boost::default_bfs_visitor {
public:
	Forest* forest;
	boost::unordered_map<Edge, SimpleEdge*>* boostToSimpleEdgeMap;

	tree_visitor(Forest* f, boost::unordered_map<Edge, SimpleEdge*>* edgeMap) : forest(f), boostToSimpleEdgeMap(edgeMap) {}

	void tree_edge(const Edge &e, const BoostGraph &g) const {
		Vertex src = source(e, g);
		Vertex tar = target(e, g);
		//CGAL::Lazy_exact_nt<CGAL::Gmpq> exactWeight = CGAL::squared_distance(g[src].pt, g[tar].pt);
		//double weight = CGAL::to_double(exactWeight); // Re-calculate weights since constraint edges have 0 weight

		// For Link(), first param is the leafmost node, 2nd param is the parent
		forest->Link(tar, src);

		SimpleEdge* correspondingEdge = (*boostToSimpleEdgeMap)[e];
		int sortedOrder = correspondingEdge->sortedOrder;
		forest->SetCost(tar, sortedOrder);

		if (SHOW_DEBUG) {
			std::cout << "Tree edge: w:" << sortedOrder << " : " << tar << " (" << g[tar].pt << ") " << src << " (" << g[src].pt << ")" << std::endl;
		}
	}
};

Forest* createLinkCutTree(BoostGraph* g, boost::unordered_map<Edge, SimpleEdge*>* boostToSimpleEdgeMap) {
	Forest* f = new Forest();
	f->Initialize(g->m_vertices.size());

	tree_visitor vis(f, boostToSimpleEdgeMap);

	// Arbitrary root vertex, s
	Vertex s = *(boost::vertices(*g).first);
	boost::breadth_first_search(*g, s, boost::visitor(vis));

	// Because weights are stored in the leafmost node, the node s represents the cost
	// for s and null node (parent of root). This is an invalid edge, so we 
	// initialize it with an invalid weight s.t. GetCost(s) is smaller than all possible costs.
	f->SetCost(s, 0);

	return f;
}

EdgeVector* checkCycles(Forest* f, EdgeVector* edgesTPrime, boost::unordered_set<SimpleEdge>* constraintSet, boost::unordered_map<SimpleEdge, SimpleEdge*>* relativeWeightSet) {

	EdgeVector* sEdges = new EdgeVector();

	for (int i = 0; i < edgesTPrime->size(); i++) {
		SimpleEdge* edgeS = (*edgesTPrime)[i];
		VertexIndex u = edgeS->u;
		VertexIndex v = edgeS->v;
		SimpleEdge* cdtEdge = (*relativeWeightSet)[(*edgeS)];
		int we = cdtEdge->sortedOrder;
		assert(we > 0);

		// f->SameEdge(u, v) determines if an edge(u, v) from TPrime already exists in Cmst. Adding an already existing edge to Cmst 
		// is guaranteed to not form a cycle and thus is skipped
		if (f->SameEdge(u, v)) {
			continue;
		}

		// f->SameTree(u, v) determines if an edge(u, v) from TPrime belong to the same connected component in Cmst.
		if (f->SameTree(u, v)) {
			Forest::NodeId lca = f->LCA(u, v);

			// LCA weight needs to be ignored (zero weight) because all weights are stored in the leafmost 
			// node, i.e. weight on node lca is the edge weight for edge (lca, parent lca)
			int lcaWeight = f->GetCost(lca);
			f->SetCost(lca, 0);

			bool hasParent = f->FindParent(lca) != Forest::UNDEFINED_NODEID;
			Forest::NodeId parentOfLca = Forest::UNDEFINED_NODEID;
			if (hasParent) {
				parentOfLca = f->Cut(lca);
			}

			Forest::NodeId max_u = f->FindMax(u);
			Forest::Cost cost_u = f->GetCost(max_u);
			while (we <= cost_u) {
				f->SetCost(max_u, 0);

				// Handling non-unique costs in the algorithm requires extra work. If there is a cycle C
				// where every edge has equal weight, removing any edge creates an MST. Assume there is an edge e in C that is a constraint.
				// There are a few possible outcomes:
				// 1) The MST algo picks e for TPrime (then e will not exist in S)
				// 2) The MST algo does not pick e for TPrime (then e will still not be added to S since it has equal weight as all other edges)
				// Our solution to this is to check all cycles when they have equal or lesser weight (instead of lesser weight exclusively) and add
				// an edge to S if it originally was a constrained edge
				if (constraintSet->count(SimpleEdge(max_u, f->FindParent(max_u), 0)) > 0) {
					SimpleEdge* edgeS = new SimpleEdge(max_u, f->FindParent(max_u), 0);
					sEdges->push_back(edgeS);
				}

				max_u = f->FindMax(u);
				cost_u = f->GetCost(max_u);
			}

			Forest::NodeId max_v = f->FindMax(v);
			Forest::Cost cost_v = f->GetCost(max_v);
			while (we <= cost_v) {
				f->SetCost(max_v, 0);

				if (constraintSet->count(SimpleEdge(max_v, f->FindParent(max_v), 0)) > 0) {
					SimpleEdge* edgeS = new SimpleEdge(max_v, f->FindParent(max_v), 0);
					sEdges->push_back(edgeS);
				}

				max_v = f->FindMax(v);
				cost_v = f->GetCost(max_v);
			}

			if (hasParent) {
				f->Link(lca, parentOfLca);
			}
			f->SetCost(lca, lcaWeight);
		}
	}

	return sEdges;
}

void computeCmst(
	VertexVector* verticesF,
	EdgeVector* edgesF,
	EdgeVector** NewEdgesF,
	EdgeVector** TPrime,
	EdgeVector** Cmst,
	EdgeVector** S)
{
	// Input: plane forest F = (V, E)
	// Output: minimum set S ⊆ E of constraints such that newF ⊆ CMST(V, S)
	// if (SHOW_DEBUG) { printGraph("Input graph F = (V, E)", F, true); }

	boost::chrono::high_resolution_clock::time_point start;
	boost::chrono::high_resolution_clock::time_point end;
	boost::chrono::milliseconds duration(0);
	boost::chrono::milliseconds total(0);

	// Compute CDT(F)
	start = boost::chrono::high_resolution_clock::now();
	CDT* cdt = computeCdt(verticesF, edgesF);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("CDT(F)", duration);
	if (SHOW_DEBUG) { printCdtInfo(cdt); }

	// Replace F with NewF
	start = boost::chrono::high_resolution_clock::now();
	(*NewEdgesF) = newConstraintSetFromCdt(cdt, verticesF);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("F -> NewF", duration);
	//if (SHOW_DEBUG) { printGraph("F -> NewF by splitting collinear constraints", *NewF); }

	// CDT -> EdgeVector
	start = boost::chrono::high_resolution_clock::now();
	EdgeVector* cdtEdgeVector = convertCdtToGraph(verticesF, cdt);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("CDT -> EdgeVector", duration);
	//if (SHOW_DEBUG) { printGraph("F -> NewF by splitting collinear constraints", *NewF); }

	// Recompute sort order (relative weight)
	start = boost::chrono::high_resolution_clock::now();
	boost::unordered_map<SimpleEdge, SimpleEdge*>* relativeWeightSet = new boost::unordered_map<SimpleEdge, SimpleEdge*>();
	cdtEdgeVector = sortByWeight(verticesF, cdtEdgeVector);
	// Assign sorted order
	for (int i = 0; i < cdtEdgeVector->size(); i++) {
		SimpleEdge* e = (*cdtEdgeVector)[i];
		e->sortedOrder = i + 1;
		(*relativeWeightSet)[(*e)] = e;
	}
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Recompute sort order", duration);

	// Compute T' = MST(CDT(NewF)) (Best case MST)
	start = boost::chrono::high_resolution_clock::now();
	(*TPrime) = computeCustomMst(verticesF, cdtEdgeVector);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("T' = MST(CDT(F))", duration);
	//if (SHOW_DEBUG) { printGraph("T' = MST(CDT(F)) Best case MST", *TPrime); }

	// Compute CMST(NewF) = MST(CDT◦(NewF))
	start = boost::chrono::high_resolution_clock::now();
	boost::unordered_set<SimpleEdge>* contraintEdgeSet = createSimpleEdgeSet(*NewEdgesF);
	(*Cmst) = computeCustomMst(verticesF, cdtEdgeVector, contraintEdgeSet);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("CMST(F) = MST(CDT◦(F))", duration);
	//if (SHOW_DEBUG) { printGraph("CMST(F) = MST(CDT◦(F))", *Cmst); }

	// Generate Boost graph for CMST
	start = boost::chrono::high_resolution_clock::now();
	boost::unordered_map<Edge, SimpleEdge*>* boostToSimpleEdgeMap;
	BoostGraph* CmstBg = convertVectorToGraph(verticesF, (*Cmst), &boostToSimpleEdgeMap, relativeWeightSet);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Generate Boost graph for CMST", duration);

	// Create DynamicTree for CMST(NewF)
	start = boost::chrono::high_resolution_clock::now();
	Forest* dt = createLinkCutTree(CmstBg, boostToSimpleEdgeMap);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Dynamic Tree construction", duration);

	// Check for cycles and construct constraint set s
	start = boost::chrono::high_resolution_clock::now();
	*S = checkCycles(dt, (*TPrime), contraintEdgeSet, relativeWeightSet);
	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	total += duration;
	printDuration("Check cycles", duration);
	//if (SHOW_DEBUG) { printGraph("Contraint set S", *S); }

	printDuration("Total", total);

	delete dt;
	//delete CmstBg;
	delete contraintEdgeSet;
	deleteEdgeVector(cdtEdgeVector);
	delete relativeWeightSet;
	delete cdt;
}

bool containsEdge(boost::unordered_set<SimpleEdge>* edgeSet, SimpleEdge* edge) {
	SimpleEdge se(edge->u, edge->v, 0);
	return edgeSet->count(se) > 0;
}

// True if A a subgraph of B
bool isSubgraph(VertexVector* vertices, EdgeVector* a, EdgeVector* b) {
	boost::unordered_set<SimpleEdge>* bEdgeSet = createSimpleEdgeSet(b);

	// Iterate through the edges
	for (int i = 0; i < a->size(); i++) {
		SimpleEdge* edge = (*a)[i];
		if (!containsEdge(bEdgeSet, edge)) {
			CGALPoint* u = (*vertices)[edge->u];
			CGALPoint* v = (*vertices)[edge->v];
			CGAL::Lazy_exact_nt<CGAL::Gmpq> exactWeight = CGAL::squared_distance(*u, *v);
			std::cout << "b contains edge not in a: " << CGAL::to_double(exactWeight) << " (" << *u << ") (" << *v << ")" << std::endl;
			return false;
		}
	}

	return true;
}

EdgeVector* graphOmitEdge(EdgeVector* edgesS, int omitIndex) {
	EdgeVector* newS = new EdgeVector();
	newS->reserve(edgesS->size());

	int index = 0;
	for (int i = 0; i < edgesS->size(); i++) {
		if (index++ == omitIndex) {
			continue;
		}

		SimpleEdge* edge = (*edgesS)[i];
		SimpleEdge* newEdge = new SimpleEdge(edge->u, edge->v, 0);
		newS->push_back(newEdge);
	}

	return newS;
}

bool isCdtSubgraph(VertexVector* vertices, EdgeVector* edgesF, EdgeVector* edgesS) {
	CDT* cdtS = computeCdt(vertices, edgesS); // CDT of mimimum edge constraint
	EdgeVector* ev_cdtS = convertCdtToGraph(vertices, cdtS);
	bool res = isSubgraph(vertices, edgesF, ev_cdtS);

	delete ev_cdtS;
	delete cdtS;

	return res;
}

bool isCmstSubgraph(VertexVector* vertices, EdgeVector* edgesF, EdgeVector* edgesS) {
	CDT* cdtS = computeCdt(vertices, edgesS); // CDT of mimimum edge constraint
	EdgeVector* ev_cdtS = convertCdtToGraph(vertices, cdtS);
	EdgeVector* cmstS = computeCustomMst(vertices, ev_cdtS, createSimpleEdgeSet(edgesS));
	bool res = isSubgraph(vertices, edgesF, cmstS);

	delete cmstS;
	delete ev_cdtS;
	delete cdtS;

	return res;
}

bool isMinimal(VertexVector* vertices, EdgeVector* edgesF, EdgeVector* edgesS) {
	// Removal of any edge of S should result in F !⊆ CMST(V, S)
	for (int i = 0; i < edgesS->size(); i++) {
		EdgeVector* omittedS = graphOmitEdge(edgesS, i);
		bool sub = isCmstSubgraph(vertices, edgesF, omittedS);
		delete omittedS;

		if (sub) {
			return false;
		}
	}
	return true;
}

bool isContainedIn(EdgeVector* edgesF, EdgeVector* edgesS, EdgeVector* edgesTPrime) {
	boost::unordered_set<SimpleEdge>* sEdgeSet = createSimpleEdgeSet(edgesS);
	boost::unordered_set<SimpleEdge>* tPrimeEdgeSet = createSimpleEdgeSet(edgesTPrime);

	for (int i = 0; i < edgesF->size(); i++) {
		SimpleEdge* edge = (*edgesF)[i];
		if (containsEdge(sEdgeSet, edge))
			continue;
		else if (containsEdge(tPrimeEdgeSet, edge))
			continue;
		else
			return false;
	}

	return true;
}

// TODO:
// Move the PERFORM_RESTRICTION_CHECKS into it's own project
int main(int argc, char* argv[]) {
	const char* vertFile = (argc > 2) ? argv[1] : NULL;
	const char* edgeFile = (argc > 2) ? argv[2] : NULL;

	VertexVector* vertices = NULL;
	EdgeVector* edges = NULL;

	if (vertFile == NULL || edgeFile == NULL) {
		// Random graph
		//createRandomPlaneForest(1000, 1000, 100, &vertices, &edges);
		//createRandomPlaneForest(50, 1000, 25, &vertices, &edges);
		//createRandomNearTriangulation(1000, 1000, &vertices, &edges);
		createRandomNearTriangulation(100, 1000, &vertices, &edges);
		//printGDFGraph("D:\\g\\results\\graph examples\\randPlaneInDisc_1000_1000_mst.gdf", vertices, edges);
	}
	else {
		// Load graph from file
		// E.g. D:\g\data\HI.nodes D:\g\data\HI.edges
		parseGraph(vertFile, edgeFile, &vertices, &edges);
		//printGDFGraph(vertices, edges);
	}

	/*vertices = new VertexVector();
	vertices->push_back(new CGALPoint(0, 0));
	vertices->push_back(new CGALPoint(1, 0));
	vertices->push_back(new CGALPoint(1, 1));

	edges = new EdgeVector();
	edges->push_back(new SimpleEdge(0, 1, 0));
	edges->push_back(new SimpleEdge(2, 0, 0));*/

	EdgeVector* NewEdges = NULL;
	EdgeVector* TPrime = NULL;
	EdgeVector* Cmst = NULL;
	EdgeVector* S = NULL;

	computeCmst(vertices, edges, &NewEdges, &TPrime, &Cmst, &S);

	//printGDFGraph("D:\\g\\results\\graph examples\\in.gdf", vertices, NewEdges);
	//printGDFGraph("D:\\g\\results\\graph examples\\out.gdf", vertices, S);

	// Edge Ratio
	std::cout << "Edges in E: " << NewEdges->size() << " Edges in S: " << S->size() << " Ratio: " << (double)((double)S->size() / (double)NewEdges->size()) << std::endl;

	// Validatation:
	// S ⊆ E
	if (!isSubgraph(vertices, S, NewEdges)) {
		std::cout << "Error: isSubgraph is false" << std::endl;
	}

	// For each e in NewF, if e not in S, then e in T'
	if (!isContainedIn(NewEdges, S, TPrime)) {
		std::cout << "Error: isContainedIn is false" << std::endl;
	}

	// F ⊆ CDT(V, S) (Note: CMST is contained in CDT)
	if (!isCdtSubgraph(vertices, NewEdges, S)) {
		std::cout << "Error: isCdtSubgraph is false" << std::endl;
	}

	// S ⊆ E s.t. F ⊆ CMST(V, S). Note: CMST(G) = MST(CVG(G)) = MST(CDT◦(G)) (where CDT◦(G) = CDT(G) when all edges in G have 0 weight)
	// Notice: CMST(V, S) is a spanning graph, that is there is at least 1 edge that connects every vertex. So F, the constraint set should be a subset of CMST(V, S)
	// In other words, we want to find the smallest subset S of edges of F such that
	// CMST(F) is equal to CMST(V, S), although the weights of the two trees may
	// be different.
	// NewF ⊆ CMST(V, S)
	if (!isCmstSubgraph(vertices, NewEdges, S)) {
		std::cout << "Error: isCmstSubgraph is false" << std::endl;
	}

	// Removal of any edge of S should result in F !⊆ CMST(V, S) (WARNING: very very slow)
	/*if (!isMinimal(vertices, NewEdges, S)) {
		std::cout << "Error: isMinimal is false" << std::endl;
		}*/

	deleteEdgeVector(S);
	deleteEdgeVector(Cmst);
	deleteEdgeVector(TPrime);
	deleteEdgeVector(NewEdges);
	deleteEdgeVector(edges);
	deleteVerticesVector(vertices);

	return 0;
}