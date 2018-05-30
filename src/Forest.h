// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com> and Brown University
// Mofified by eduong
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef FOREST_H_
#define FOREST_H_

#include "dtree/tree.h"
#include "dtree/tree-inl.h"

class Forest {
public:
	// don't change Cost to an unsigned type; that breaks the Min aggregate
	// read the appendix before changing Cost to a floating-point type
	typedef int Cost;
	const Cost UNINITIALIZED = -1;

	typedef size_t NodeId;
	const static size_t UNDEFINED_NODEID = -1;

private:
	typedef dtree::WithAncAggr < dtree::Max<Cost>,
		dtree::WithAncValue < dtree::Add<Cost>,
		dtree::Begin<> > > WithMaxCost;
	typedef dtree::EndTree<WithMaxCost> Node;

public:
	Forest() : node_(NULL) {}
	~Forest()  {
		delete[] node_;
	}

	// Creates a forest with node ids 0..size-1
	void Initialize(size_t size)  {
		Node* old_node = node_;

		// avoids a double delete[] if new[] throws
		node_ = NULL;

		delete[] old_node;

		node_size = size;
		node_ = new Node[node_size];
		for (NodeId i = 0; i < node_size; i++) {
			SetCost(i, 0);
		}
	}

	// On simulating edge values http://www.davideisenstat.com/dtree/#simulating-edge-values
	// Node u is child of v, hence it will store the cost of the edge and the EdgeId map at
	// index u contains an edge identifier
	Cost GetCost(NodeId u) {
		return WithMaxCost::Value(&node_[u]);
	}

	void SetCost(NodeId u, Cost x) {
		WithMaxCost::SetValue(&node_[u], x);
	}

	NodeId LCA(NodeId u, NodeId v) {
		return dtree::LeafmostCommonAnc(&node_[u], &node_[v]) - node_;
	}

	NodeId FindRoot(NodeId u) {
		return Root(&node_[u]) - node_;
	}

	NodeId FindParent(NodeId u) {
		Node* parent = dtree::Parent(&node_[u]);
		if (parent) {
			return parent - node_;
		}
		return UNDEFINED_NODEID;
	}

	NodeId FindMax(NodeId u) {
		// There seems to be a bug where this function often returns an invalid address when the tree is large.
		Node* rootmostAnc = WithMaxCost::FindRootmostAnc(&node_[u], dtree::GreaterEqual(WithMaxCost::AggrAnc(&node_[u])));
		return rootmostAnc - node_;
	}

	NodeId FindMax2(NodeId u) {
		NodeId root = FindRoot(u);
		NodeId n = u;
		double maxCost = 0;
		NodeId maxNode = u;

		while (n != root) {
			Cost nCost = GetCost(n);
			if (nCost > maxCost) {
				maxCost = nCost;
				maxNode = n;
			}
			n = FindParent(n);
		}

		return maxNode;
	}

	void Link(NodeId u, NodeId v)  {
		if (u == UNDEFINED_NODEID || v == UNDEFINED_NODEID) {
			return;
		}
		dtree::Link(&node_[u], &node_[v]);
	}

	NodeId Cut(NodeId u)  {
		Node* parent = NULL;
		if (parent = dtree::Parent(&node_[u])) { // Cut only if there is a parent to cut
			// "If u is nonnull and has a parent p, removes the edge from u to p and returns an ancestor of p".
			// We would actually like to return p and not p's ancestor
			//return dtree::Cut(&node_[u]) - node_;
			dtree::Cut(&node_[u]);
			return parent - node_;
		}
		return UNDEFINED_NODEID;
	}

	bool SameTree(NodeId u, NodeId v) {
		return dtree::SameTree(&node_[u], &node_[v]);
	}

	bool SameEdge(NodeId u, NodeId v) {
		return u == (dtree::Parent(&node_[v]) - node_) || v == (dtree::Parent(&node_[u]) - node_);
	}

private:
	Node* node_;
	size_t node_size;

	// disallows the copy constructor and the assignment operator
	Forest(const Forest&);
	void operator=(const Forest&);
};

#endif  // FOREST_H_