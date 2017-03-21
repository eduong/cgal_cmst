// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef DTREE_SEQ_INL_H_
#define DTREE_SEQ_INL_H_

namespace dtree {

static const int kLeft = kLw;
static const int kRight = kRw;

template<typename Node>
struct SolidPrimitives {
  static void AddValues(Node* node, const Node* parent) {
    node->SolidAddValues(parent);
  }
  static void SubtractValues(Node* node, const Node* parent) {
    node->SolidSubtractValues(parent);
  }
  static void ResetAggrs(Node* node) { node->SolidResetAggrs(); }
  static void MergeAggrs(Node* node, const Node* solid) {
    node->SolidMergeSolidAggrs(solid);
  }
  static int Flipped(const Node* node) { return node->SolidFlipped(); }
  static Node* child(const Node* node, int i) { return node->solid(i); }
  static void set_child(Node* node, int i, Node* child_i) {
    node->set_solid(i, child_i);
  }
  static bool Index(const Node* node, const Node* parent, int* i) {
    return node->SolidIndex(parent, i);
  }
};

template<typename Node, typename Prim>
struct SplayFunctions {
  static void RotateUp(Node* u, int index_of_e, Node* p, int index_of_u) {
    //    p     u    .
    //   /       \   .
    //  u   ==>   p  .
    //   \       /   .
    //    e     e    .
    Node* e = Prim::child(u, index_of_e);
    Prim::set_child(u, index_of_e, p);
    Prim::SubtractValues(p, u);
    Prim::ResetAggrs(p);
    Prim::set_child(p, index_of_u, e);
    if (e) {
      Prim::SubtractValues(e, p);
      e->set_parent(p);
      Prim::MergeAggrs(p, e);
    }
    if (Node* w = Prim::child(p, index_of_u ^ 1)) Prim::MergeAggrs(p, w);
  }
  template<bool finishing_search>
  static Node* Splay(Node* u);
  static Node* Dir(int dir, Node* u);
};

template<typename Node, typename Prim>
template<bool finishing_search>
Node* SplayFunctions<Node, Prim>::Splay(Node* u) {
  Node* p = u->parent();
  int index_of_u;
  if (!Prim::Index(u, p, &index_of_u)) return u;
  Node* gp = p->parent();
  int index_of_p;
  if (!Prim::Index(p, gp, &index_of_p)) {
    //      p        u    .
    //     /        / \   .
    //    u   ==>  c   p  .
    //   / \          /   .
    //  c   e        e    .
    if (!finishing_search) Prim::AddValues(u, p);
    int index_of_c = index_of_u ^ Prim::Flipped(p) ^ Prim::Flipped(u);
    RotateUp(u, index_of_c ^ 1, p, index_of_u);
    u->set_parent(gp);
    p->set_parent(u);
    return p;
  }
  Node* old_root;
  while (true) {
    if (!finishing_search) {
      Prim::AddValues(p, gp);
      Prim::AddValues(u, p);
    }
    int index_of_c = index_of_u ^ Prim::Flipped(p) ^ Prim::Flipped(u);
    u->set_parent(gp->parent());
    if ((index_of_p ^ Prim::Flipped(gp)) == (index_of_u ^ Prim::Flipped(p))) {
      //        gp        u       .
      //       /         / \      .
      //      p         c   p     .
      //     / \   ==>     / \    .
      //    u   w         e   gp  .
      //   / \               /    .
      //  c   e             w     .
      RotateUp(p, index_of_u ^ 1, gp, index_of_p);
      gp->set_parent(p);
      RotateUp(u, index_of_c ^ 1, p, index_of_u);
    } else {
      //    gp                   .
      //     \            u      .
      //      p         /   \    .
      //     /   ==>  gp      p  .
      //    u          \     /   .
      //   / \          c   e    .
      //  c   e                  .
      RotateUp(u, index_of_c ^ 1, p, index_of_u);
      RotateUp(u, index_of_c, gp, index_of_p);
    }
    p = u->parent();
    if (!Prim::Index(gp, p, &index_of_u)) {
      old_root = gp;
      break;
    }
    gp = p->parent();
    if (Prim::Index(p, gp, &index_of_p)) continue;
    if (!finishing_search) Prim::AddValues(u, p);
    index_of_c = index_of_u ^ Prim::Flipped(p) ^ Prim::Flipped(u);
    RotateUp(u, index_of_c ^ 1, p, index_of_u);
    u->set_parent(gp);
    old_root = p;
    break;
  }
  if (Node* c = Prim::child(u, 0)) c->set_parent(u);
  if (Node* e = Prim::child(u, 1)) e->set_parent(u);
  return old_root;
}

template<typename Node, typename Prim>
Node* SplayFunctions<Node, Prim>::Dir(int dir, Node* u) {
  int index_of_c = dir ^ Prim::Flipped(u);
  while (Node* c = Prim::child(u, index_of_c)) {
    Prim::AddValues(c, u);
    int index_of_gc = dir ^ Prim::Flipped(c);
    Node* gc = Prim::child(c, index_of_gc);
    if (!gc) {
      //    u       c    .
      //   /         \   .
      //  c    ==>    u  .
      //   \         /   .
      //    ge      ge   .
      RotateUp(c, index_of_gc ^ 1, u, index_of_c);
      u->set_parent(c);
      u = c;
      break;
    }
    //         u         gc   .
    //        /         / \   .
    //       c        ggc  u  .
    //      /    ==>      /   .
    //     gc            c    .
    //    / \           /     .
    //  ggc gge       gge     .
    Prim::AddValues(gc, c);
    int index_of_ggc = dir ^ Prim::Flipped(gc);
    RotateUp(gc, index_of_ggc ^ 1, c, index_of_gc);
    RotateUp(gc, index_of_ggc ^ 1, u, index_of_c);
    u->set_parent(gc);
    u = gc;
    index_of_c = index_of_ggc;
    c = Prim::child(u, index_of_c);
  }
  return u;
}

template<typename Base>
template<typename Forwarder>
inline void EndSeq<Base>::CopyFrom(const EndSeq& node, Forwarder forward) {
  Base::operator=(node);
  set_parent(node.parent() ? forward(this, &node, node.parent()) : NULL);
  set_solid(0, node.solid(0) ? forward(this, &node, node.solid(0)) : NULL);
  set_solid(1, node.solid(1) ? forward(this, &node, node.solid(1)) : NULL);
}

template<typename Base>
inline EndSeq<Base>* EndSeq<Base>::Expose() {
  typedef SplayFunctions<EndSeq, SolidPrimitives<EndSeq> > Solid;
  Solid::template Splay<false>(this);
  return this;
}

template<typename Base>
inline void EndSeq<Base>::FinishSearch() {
  typedef SplayFunctions<EndSeq, SolidPrimitives<EndSeq> > Solid;
  Solid::template Splay<true>(this);
}

template<typename Node>
inline void UpdateAggrs(Node* u) {
  u->SolidResetAggrs();
  if (Node* c = u->solid(0)) u->SolidMergeSolidAggrs(c);
  if (Node* e = u->solid(1)) u->SolidMergeSolidAggrs(e);
}

template<typename Semigroup, typename Self, typename Node, typename Predicate>
class Search {
 public:
  Search(Node* node, const Predicate& predicate)
      : node_(node->Expose()),
        predicate_(predicate),
        cum_aggr_(Semigroup::empty_aggr()) {
    assert(!predicate(Semigroup::empty_aggr()));
  }
  ~Search() { node()->FinishSearch(); }
  Node* FindDirmostSeq(int dir) {
    while (true) {
      int i = dir ^ node()->SolidFlipped();
      if (TestSolid(i)) continue;
      if (TestNode()) return node();
      if (TestSolid(i ^ 1)) continue;
      return NULL;
    }
  }

 protected:
  Node* node() const { return node_; }
  void set_node(Node* node) { node_ = node; }
  bool Test(typename Semigroup::Type aggr) {
    typename Semigroup::Type new_cum_aggr =
        Semigroup::CombineAggrs(cum_aggr_, aggr);
    if (predicate_(new_cum_aggr)) return true;
    cum_aggr_ = new_cum_aggr;
    return false;
  }
  bool TestAggr(Node* child) { return Test(child->Self::aggr()); }
  bool TestSolid(int i) {
    Node* solid = node()->solid(i);
    if (!solid) return false;
    solid->SolidAddValues(node());
    if (!TestAggr(solid)) {
      solid->SolidSubtractValues(node());
      return false;
    }
    set_node(solid);
    return true;
  }
  bool TestNode() { return Test(node()->Self::singleton_aggr()); }

 private:
  Node* node_;
  const Predicate& predicate_;
  typename Semigroup::Type cum_aggr_;
  Search(const Search&);
  void operator=(const Search&);
};

////////////////////////////////////////

template<typename Node>
inline Node* Dir(int dir, Node* node) {
  typedef SplayFunctions<Node, SolidPrimitives<Node> > Solid;
  if (!node) return NULL;
  node = Solid::Dir(dir, node->Expose());
  node->set_parent(NULL);
  return node;
}

template<typename Node>
inline Node* Left(Node* node) {
  return Dir(kLeft, node);
}

template<typename Node>
inline Node* Right(Node* node) {
  return Dir(kRight, node);
}

template<typename Node>
inline bool Connected(Node* node1, Node* node2) {
  if (node1 == node2) return true;
  if (!(node1 && node2)) return false;
  node1->Expose();
  node2->Expose();
  return node1->parent();
}

template<typename Node>
inline bool SameSeq(Node* node1, Node* node2) {
  return Connected(node1, node2);
}

template<typename Node>
inline Node* CutDirOf(int dir, Node* node2) {
  if (!node2) return NULL;
  int i = dir ^ node2->Expose()->SolidFlipped();
  Node* node1 = node2->solid(i);
  if (!node1) return NULL;
  node2->set_solid(i, NULL);
  node1->SolidAddValues(node2);
  node1->set_parent(NULL);
  return node1;
}

template<typename Node>
inline Node* CutLeftOf(Node* node2) {
  return CutDirOf(kLeft, node2);
}

template<typename Node>
inline Node* CutRightOf(Node* node2) {
  return CutDirOf(kRight, node2);
}

template<typename Node>
inline Node* LinkDirOf(int dir, Node* node1, Node* node2) {
  if (!node1) return node2;
  if (!node2) return node1;
  UpdateAggrs(node1->Expose());
  node2->Expose();
  assert(!node1->parent());
  node2 = Dir(dir, node2);
  node2->set_solid(dir ^ node2->SolidFlipped(), node1);
  node1->SolidSubtractValues(node2);
  node1->set_parent(node2);
  return node2;
}

template<typename Node>
inline Node* LinkLeftOf(Node* node1, Node* node2) {
  return LinkDirOf(kLeft, node1, node2);
}

template<typename Node>
inline Node* LinkRightOf(Node* node1, Node* node2) {
  return LinkDirOf(kRight, node1, node2);
}

template<typename Node>
class ScopedCutDirOf {
 public:
  ScopedCutDirOf(int dir, Node* node2)
      : dir_(dir), excl_part_(CutDirOf(dir, node2)), incl_part_(node2) {
  }
  ~ScopedCutDirOf() { LinkDirOf(dir(), excl_part(), incl_part()); }
  int dir() const { return dir_; }
  Node* excl_part() const { return excl_part_; }
  Node* incl_part() const { return incl_part_; }
 private:
  int dir_;
  Node* excl_part_;
  Node* incl_part_;
  ScopedCutDirOf(const ScopedCutDirOf&);
  void operator=(const ScopedCutDirOf&);
};

template<typename Node>
class ScopedCutLeftOf : public ScopedCutDirOf<Node> {
 public:
  explicit ScopedCutLeftOf(Node* node2) : ScopedCutDirOf<Node>(kLeft, node2) {}
 private:
  ScopedCutLeftOf(const ScopedCutLeftOf&);
  void operator=(const ScopedCutLeftOf&);
};

template<typename Node>
class ScopedCutRightOf : public ScopedCutDirOf<Node> {
 public:
  explicit ScopedCutRightOf(Node* node2)
      : ScopedCutDirOf<Node>(kRight, node2) {
  }
 private:
  ScopedCutRightOf(const ScopedCutRightOf&);
  void operator=(const ScopedCutRightOf&);
};

template<typename Node>
inline Node* Dirward(int dir, Node* node) {
  ScopedCutDirOf<Node> dir_cut(dir, node);
  return Dir(dir ^ 1, dir_cut.excl_part());
}

template<typename Node>
inline Node* Leftward(Node* node) {
  return Dirward(kLeft, node);
}

template<typename Node>
inline Node* Rightward(Node* node) {
  return Dirward(kRight, node);
}

template<typename Group, typename Base>
template<typename Node>
inline typename Group::Type WithValue<Group, Base>::Value(Node* node) {
  typedef WithValue Self;
  assert(node);
  return node->Expose()->Self::value();
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SetValue(Node* node,
                                             typename Group::Type value) {
  typedef WithValue Self;
  assert(node);
  ScopedCutLeftOf<Node> left_cut(node);
  ScopedCutRightOf<Node> right_cut(node);
  node->Self::set_value(value);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::AddToSeq(Node* node,
                                             typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  node->Expose()->Self::add_to_value(delta);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SubtractFromSeq(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  node->Expose()->Self::subtract_from_value(delta);
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type WithAggr<Semigroup, Base>::AggrSeq(
    Node* node) {
  typedef WithAggr Self;
  if (!node) return Semigroup::empty_aggr();
  UpdateAggrs(node->Expose());
  return node->Self::aggr();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithAggr<Semigroup, Base>::FindDirmostSeq(int dir,
                                                Node* node,
                                                const Predicate& predicate) {
  if (!node) return NULL;
  Search<Semigroup, WithAggr, Node, Predicate> search(node, predicate);
  return search.FindDirmostSeq(dir);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindLeftmostSeq(
    Node* node,
    const Predicate& predicate) {
  return FindDirmostSeq(kLeft, node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindRightmostSeq(
    Node* node,
    const Predicate& predicate) {
  return FindDirmostSeq(kRight, node, predicate);
}

template<typename Base>
template<typename Node>
inline void WithReverseBy<Base>::ReverseSeq(Node* node) {
  typedef WithReverseBy Self;
  Self::AddToSeq(node, Group::flip_delta());
}
}  // namespace dtree
#endif  // DTREE_SEQ_INL_H_
