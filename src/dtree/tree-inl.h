// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef DTREE_TREE_INL_H_
#define DTREE_TREE_INL_H_

#include "seq-inl.h"

namespace dtree {

static const int kLeaf = kLw;
static const int kRoot = kRw;

template<bool finishing_search, typename Node>
Node* ExposeTemplate(Node* u) {
  typedef SplayFunctions<Node, SolidPrimitives<Node> > Solid;
  Node* even_d0;
  Node* even_d1;
  u->template ExposeStep<finishing_search>(&even_d0, &even_d1);
  Node* p = u->parent();
  if (!p) return NULL;
  int index_of_e = kRoot ^ u->SolidFlipped();
  Node* last_splice;
  while (true) {
    Node* odd_d0;
    Node* odd_d1;
    p->template ExposeStep<finishing_search>(&odd_d0, &odd_d1);
    int index_of_v = kLeaf ^ p->SolidFlipped();
    if (Node* v = p->solid(index_of_v)) {
      p->template DashedReplace<finishing_search>(v, even_d0, even_d1);
    } else {
      p->template DashedRemove<finishing_search>(even_d0, even_d1);
    }
    Node* gp = p->parent();
    if (!gp) {
      //    p       u    .
      //   /|        \   .
      //  v u   ==>   p  .
      //     \       /|  .
      //      e     e v  .
      if (!finishing_search) u->DottedAddValues(p);
      Solid::RotateUp(u, index_of_e, p, index_of_v);
      last_splice = p;
      break;
    }
    gp->template ExposeStep<finishing_search>(&even_d0, &even_d1);
    int index_of_q = kLeaf ^ gp->SolidFlipped();
    if (Node* q = gp->solid(index_of_q)) {
      gp->template DashedReplace<finishing_search>(q, odd_d0, odd_d1);
    } else {
      gp->template DashedRemove<finishing_search>(odd_d0, odd_d1);
    }
    Node* ggp = gp->parent();
    //    gp        u       .
    //   /|          \      .
    //  q p           p     .
    //   /|\   ==>   /|\    .
    //  v u w       e v gp  .
    //     \           /|   .
    //      e         w q   .
    if (!finishing_search) {
      p->DottedAddValues(gp);
      u->DottedAddValues(p);
    }
    Solid::RotateUp(p, index_of_v ^ 1, gp, index_of_q);
    gp->set_parent(p);
    Solid::RotateUp(u, index_of_e, p, index_of_v);
    if (!ggp) {
      last_splice = gp;
      break;
    }
    p = ggp;
  }
  p->set_parent(u);
  u->set_parent(NULL);
  return last_splice;
}

template<typename Node>
inline void HardExposeTemplate(Node* u) {
  int index_of_c = kLeaf ^ u->Expose()->SolidFlipped();
  if (Node* c = u->solid(index_of_c)) {
    u->set_solid(index_of_c, NULL);
    c->SolidToDashed(u);
    u->DashedInsert(c);
  }
}

template<typename Base>
template<typename Forwarder>
inline void EndTree<Base>::CopyFrom(const EndTree& node, Forwarder forward) {
  Base::operator=(node);
  set_parent(node.parent() ? forward(this, &node, node.parent()) : NULL);
  set_solid(0, node.solid(0) ? forward(this, &node, node.solid(0)) : NULL);
  set_solid(1, node.solid(1) ? forward(this, &node, node.solid(1)) : NULL);
}

template<typename Base>
template<bool finishing_search>
inline void EndTree<Base>::ExposeStep(EndTree** /*d0*/, EndTree** /*d1*/) {
  typedef SplayFunctions<EndTree, SolidPrimitives<EndTree> > Solid;
  Solid::template Splay<finishing_search>(this);
}

template<typename Base>
template<bool finishing_search>
inline void EndTree<Base>::DashedReplace(EndTree* c,
                                         EndTree* /*d0*/,
                                         EndTree* /*d1*/) {
  c->SolidToDashed(this);
}

template<typename Base>
template<bool finishing_search>
inline void EndTree<Base>::DashedRemove(EndTree* /*d0*/, EndTree* /*d1*/) {
}

template<typename Base>
inline void EndTree<Base>::DashedInsert(EndTree* /*c*/) {
}

template<typename Base>
inline EndTree<Base>* EndTree<Base>::Expose() {
  ExposeTemplate<false>(this);
  return this;
}

template<typename Base>
inline EndTree<Base>* EndTree<Base>::HardExpose() {
  HardExposeTemplate(this);
  return this;
}

template<typename Base>
inline void EndTree<Base>::FinishSearch() {
  typedef SplayFunctions<EndTree, SolidPrimitives<EndTree> > Solid;
  Solid::template Splay<true>(this);
}

template<typename Base>
struct EndTreeWithDesc<Base>::DottedPrimitives {
  static void AddValues(EndTreeWithDesc* node, const EndTreeWithDesc* parent) {
    node->DottedAddValues(parent);
  }
  static void SubtractValues(EndTreeWithDesc* node,
                             const EndTreeWithDesc* parent) {
    node->DottedSubtractValues(parent);
  }
  static void ResetAggrs(EndTreeWithDesc* node) { node->DottedResetAggrs(); }
  static void MergeAggrs(EndTreeWithDesc* node,
                         const EndTreeWithDesc* dotted) {
    node->DottedMergeDottedAggrs(dotted);
  }
  static int Flipped(const EndTreeWithDesc* /*node*/) { return 0; }
  static EndTreeWithDesc* child(const EndTreeWithDesc* node, int i) {
    return node->dotted(i);
  }
  static void set_child(EndTreeWithDesc* node,
                        int i,
                        EndTreeWithDesc* child_i) {
    node->set_dotted(i, child_i);
  }
  static bool Index(const EndTreeWithDesc* node,
                    const EndTreeWithDesc* parent,
                    int* i) {
    return node->DottedIndex(parent, i);
  }
};

template<typename Base>
template<typename Forwarder>
inline void EndTreeWithDesc<Base>::CopyFrom(const EndTreeWithDesc& node,
                                            Forwarder forward) {
  Base::operator=(node);
  set_parent(node.parent() ? forward(this, &node, node.parent()) : NULL);
  set_dashed(node.dashed() ? forward(this, &node, node.dashed()) : NULL);
  set_solid(0, node.solid(0) ? forward(this, &node, node.solid(0)) : NULL);
  set_solid(1, node.solid(1) ? forward(this, &node, node.solid(1)) : NULL);
  set_dotted(0, node.dotted(0) ? forward(this, &node, node.dotted(0)) : NULL);
  set_dotted(1, node.dotted(1) ? forward(this, &node, node.dotted(1)) : NULL);
}

template<typename Base>
template<bool finishing_search>
inline void EndTreeWithDesc<Base>::ExposeStep(EndTreeWithDesc** d0,
                                              EndTreeWithDesc** d1) {
  typedef EndTreeWithDesc Node;
  typedef SplayFunctions<Node, SolidPrimitives<Node> > Solid;
  EndTreeWithDesc* root = Solid::template Splay<finishing_search>(this);
  if (!parent()) return;
  if (root == this) {
    DottedSplay<finishing_search>(d0, d1);
    return;
  }
  EndTreeWithDesc* root_p = root->parent();
  root->set_parent(parent());
  if (root_p == this) {
    root->DottedAddValues(this);
    if (!finishing_search) this->DottedSubtractValues(root);
    root->DottedSplay<finishing_search>(d0, d1);
    if (!finishing_search) this->DottedAddValues(root);
    root->DottedSubtractValues(this);
  } else {
    root->DottedAddValues(root_p);
    root->DottedAddValues(this);
    if (!finishing_search) this->DottedSubtractValues(root);
    root->DottedSplay<finishing_search>(d0, d1);
    if (!finishing_search) this->DottedAddValues(root);
    root->DottedSubtractValues(this);
    root->DottedSubtractValues(root_p);
  }
  set_parent(root->parent());
  root->set_parent(root_p);
}

template<typename Base>
template<bool finishing_search>
inline void EndTreeWithDesc<Base>::DashedReplace(EndTreeWithDesc* c,
                                                 EndTreeWithDesc* d0,
                                                 EndTreeWithDesc* d1) {
  c->SolidToDashed(this);
  if (d0) {
    if (finishing_search) c->DottedAddValues(this);
    c->DottedResetAggrs();
    c->DottedLink(0, d0);
    if (d1) c->DottedLink(1, d1);
    if (finishing_search) c->DottedSubtractValues(this);
  } else if (d1) {
    if (finishing_search) c->DottedAddValues(this);
    c->DottedResetAggrs();
    c->DottedLink(1, d1);
    if (finishing_search) c->DottedSubtractValues(this);
  } else {
    c->DottedResetAggrs();
  }
  set_dashed(c);
}

template<typename Base>
template<bool finishing_search>
inline void EndTreeWithDesc<Base>::DashedRemove(EndTreeWithDesc* d0,
                                                EndTreeWithDesc* d1) {
  typedef SplayFunctions<EndTreeWithDesc, DottedPrimitives> Dotted;
  if (d1) {
    if (d0) {
      d1 = Dotted::Dir(0, d1);
      d1->DottedResetAggrs();
      d1->DottedLink(0, d0);
      if (EndTreeWithDesc* d2 = d1->dotted(1)) d1->DottedMergeDottedAggrs(d2);
    }
    set_dashed(d1);
    if (finishing_search) d1->DottedSubtractValues(this);
    d1->set_parent(this);
  } else if (d0) {
    set_dashed(d0);
    if (finishing_search) d0->DottedSubtractValues(this);
    d0->set_parent(this);
  } else {
    set_dashed(NULL);
  }
}

template<typename Base>
inline void EndTreeWithDesc<Base>::DashedInsert(EndTreeWithDesc* c) {
  c->DottedResetAggrs();
  if (EndTreeWithDesc* d = dashed()) c->DottedLink(1, d);
  set_dashed(c);
}

template<typename Base>
inline EndTreeWithDesc<Base>* EndTreeWithDesc<Base>::Expose() {
  ExposeTemplate<false>(this);
  return this;
}

template<typename Base>
inline EndTreeWithDesc<Base>* EndTreeWithDesc<Base>::HardExpose() {
  HardExposeTemplate(this);
  return this;
}

template<typename Base>
inline void EndTreeWithDesc<Base>::FinishSearch() {
  ExposeTemplate<true>(this);
}

template<typename Base>
inline EndTreeWithDesc<Base>* EndTreeWithDesc<Base>::DottedCut(int i) {
  EndTreeWithDesc* dotted_i = dotted(i);
  if (!dotted_i) return NULL;
  set_dotted(i, NULL);
  dotted_i->DottedAddValues(this);
  return dotted_i;
}

template<typename Base>
inline void EndTreeWithDesc<Base>::DottedLink(int i,
                                              EndTreeWithDesc* dotted_i) {
  set_dotted(i, dotted_i);
  dotted_i->DottedSubtractValues(this);
  this->DottedMergeDottedAggrs(dotted_i);
  dotted_i->set_parent(this);
}

template<typename Base>
template<bool finishing_search>
inline void EndTreeWithDesc<Base>::DottedSplay(EndTreeWithDesc** d0,
                                               EndTreeWithDesc** d1) {
  typedef SplayFunctions<EndTreeWithDesc, DottedPrimitives> Dotted;
  Dotted::template Splay<finishing_search>(this);
  *d0 = DottedCut(0);
  *d1 = DottedCut(1);
}

template<typename Base>
inline void EndTreeWithDesc<Base>::Splice(EndTreeWithDesc* d) {
  typedef SplayFunctions<EndTreeWithDesc, DottedPrimitives> Dotted;
  if (EndTreeWithDesc* d1 = d->DottedCut(1)) {
    if (EndTreeWithDesc* d0 = d->DottedCut(0)) {
      d1 = Dotted::Dir(0, d1);
      d1->DottedResetAggrs();
      d1->DottedLink(0, d0);
      if (EndTreeWithDesc* d2 = d1->dotted(1)) d1->DottedMergeDottedAggrs(d2);
    }
    set_dashed(d1);
    d1->set_parent(this);
  } else if (EndTreeWithDesc* d0 = d->DottedCut(0)) {
    set_dashed(d0);
    d0->set_parent(this);
  } else {
    set_dashed(NULL);
  }
}

template<typename Base>
inline void EndTreeWithDesc<Base>::AbsorbDashed(EndTreeWithDesc* v) {
  typedef SplayFunctions<EndTreeWithDesc, DottedPrimitives> Dotted;
  EndTreeWithDesc* d1 = v->dashed();
  if (!d1) return;
  v->set_dashed(NULL);
  d1->DottedAddValues(v);
  if (EndTreeWithDesc* d0 = dashed()) {
    d0->DottedAddValues(this);
    d1 = Dotted::Dir(0, d1);
    d1->DottedResetAggrs();
    d1->DottedLink(0, d0);
    if (EndTreeWithDesc* d2 = d1->dotted(1)) d1->DottedMergeDottedAggrs(d2);
  }
  set_dashed(d1);
  d1->DottedSubtractValues(this);
  d1->set_parent(this);
}

////////////////////////////////////////

template<typename Node>
inline Node* Root(Node* node) {
  return Dir(kRoot, node);
}

template<typename Node>
inline bool SameTree(Node* node1, Node* node2) {
  return Connected(node1, node2);
}

template<typename Node>
inline Node* LeafmostCommonAnc(Node* node1, Node* node2) {
  if (node1 == node2) return node1;
  if (!(node1 && node2)) return NULL;
  node1->HardExpose();
  Node* last_splice = ExposeTemplate<false>(node2);
  if (!node1->parent()) return NULL;
  return last_splice ? last_splice : node2;
}

template<typename Node>
inline Node* Cut(Node* node2) {
  return CutDirOf(kRoot, node2);
}

template<typename Node>
inline Node* Link(Node* node1, Node* node2) {
  if (!node2) return node1;
  return LinkDirOf(kLeaf, node1, node2->HardExpose());
}

template<typename Node>
inline Node* Parent(Node* node) {
  return Dirward(kRoot, node);
}

template<typename Node>
Node* Leaf(Node* node) {
  typedef SplayFunctions<Node, SolidPrimitives<Node> > Solid;
  if (!node) return NULL;
  Node* u = Solid::Dir(kLeaf, node->Expose());
  while (Node* d = u->dashed()) {
    u->Splice(d);
    d->DottedAddValues(u);
    d = Solid::Dir(kLeaf, d);
    Node* gd = d->dashed();
    if (!gd) {
      //  u         d    .
      //  |          \   .
      //  d    ==>    u  .
      //   \         /   .
      //    ge      ge   .
      Solid::RotateUp(d, kRoot ^ d->SolidFlipped(),
                      u, kLeaf ^ u->SolidFlipped());
      u->set_parent(d);
      d->set_parent(NULL);
      return d;
    }
    d->Splice(gd);
    gd->DottedAddValues(d);
    gd = Solid::Dir(kLeaf, gd);
    int i = kRoot ^ gd->SolidFlipped();
    //  u            gd   .
    //  |             \   .
    //  d              u  .
    //  |     ==>     /   .
    //  gd           d    .
    //   \          /     .
    //   gge      gge     .
    Solid::RotateUp(gd, i, d, kLeaf ^ d->SolidFlipped());
    Solid::RotateUp(gd, i, u, kLeaf ^ u->SolidFlipped());
    u->set_parent(gd);
    u = gd;
  }
  u->set_parent(NULL);
  return u;
}

template<typename Node>
inline Node* Child(Node* node) {
  if (!node) return NULL;
  int index_of_c = kLeaf ^ node->Expose()->SolidFlipped();
  if (!node->solid(index_of_c)) {
    //  node       node  .
    //   |   ==>   /     .
    //   d        d      .
    Node* d = node->dashed();
    if (!d) return NULL;
    node->Splice(d);
    d->DottedAddValues(node);
    node->set_solid(index_of_c, d);
    d->SolidSubtractValues(node);
  }
  return Dirward(kLeaf, node);
}

template<typename Node>
Node* ContractDirward(int dir, Node* node) {
  Node* a = CutDirOf(kRoot, node);
  if (!a) return NULL;
  Node* p = Dir(kLeaf, a);
  Node* u;
  Node* v;
  if (dir == kLeaf) {
    u = node;
    v = p;
  } else if (dir == kRoot) {
    u = p;
    v = node;
  } else {
    assert(false);
  }
  LinkDirOf(dir ^ 1, CutDirOf(dir ^ 1, v), u);
  u->AbsorbDashed(v);
  return p;
}

template<typename Node>
inline Node* ContractLeafward(Node* node) {
  return ContractDirward(kLeaf, node);
}

template<typename Node>
inline Node* ContractRootward(Node* node) {
  return ContractDirward(kRoot, node);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::AddToAnc(Node* node,
                                             typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  node->Self::add_to_value(Group::MinusFilter(delta, delta));
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SubtractFromAnc(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  node->Self::subtract_from_value(Group::MinusFilter(delta, delta));
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::AddToProperAnc(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  if (root_cut.excl_part()) {
    root_cut.excl_part()->Self::add_to_value(Group::MinusFilter(delta, delta));
  }
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SubtractFromProperAnc(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  if (root_cut.excl_part()) {
    root_cut.excl_part()->Self::subtract_from_value(Group::MinusFilter(delta,
                                                                       delta));
  }
}

template<typename Node>
class ScopedCutDashed {
 public:
  explicit ScopedCutDashed(Node* node2)
      : excl_part_(node2->dashed()), incl_part_(node2) {
    if (excl_part()) excl_part()->DottedAddValues(node2);
  }
  ~ScopedCutDashed() {
    if (excl_part()) excl_part()->DottedSubtractValues(incl_part());
  }
  Node* excl_part() const { return excl_part_; }
  Node* incl_part() const { return incl_part_; }
 private:
  Node* excl_part_;
  Node* incl_part_;
  ScopedCutDashed(const ScopedCutDashed&);
  void operator=(const ScopedCutDashed&);
};

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SetAncDescValue(
    Node* node,
    typename Group::Type value) {
  typedef WithValue Self;
  assert(node);
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  ScopedCutDashed<Node> dashed_cut(node);
  node->Self::set_value(value);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::AddToTree(Node* node,
                                              typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  node->Expose()->Self::add_to_value_filter(delta);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SubtractFromTree(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  node->Expose()->Self::subtract_from_value_filter(delta);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::AddToDesc(Node* node,
                                              typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  node->Self::add_to_value_filter(delta);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SubtractFromDesc(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  node->Self::subtract_from_value_filter(delta);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::AddToProperDesc(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  {
    ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
    if (leaf_cut.excl_part()) {
      leaf_cut.excl_part()->Self::add_to_value_filter(delta);
    }
  }
  ScopedCutDashed<Node> dashed_cut(node);
  if (dashed_cut.excl_part()) {
    dashed_cut.excl_part()->Self::add_to_value_filter(delta);
  }
}

template<typename Group, typename Base>
template<typename Node>
inline void WithValue<Group, Base>::SubtractFromProperDesc(
    Node* node,
    typename Group::Type delta) {
  typedef WithValue Self;
  if (!node) return;
  {
    ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
    if (leaf_cut.excl_part()) {
      leaf_cut.excl_part()->Self::subtract_from_value_filter(delta);
    }
  }
  ScopedCutDashed<Node> dashed_cut(node);
  if (dashed_cut.excl_part()) {
    dashed_cut.excl_part()->Self::subtract_from_value_filter(delta);
  }
}

template<typename Type, typename Base>
template<typename Node>
inline Type WithStaticValue<Type, Base>::Value(Node* node) {
  typedef WithValue<Nop<Type>, Base> Base1;
  assert(node);
  return node->Base1::value();
}

template<typename Type, typename Base>
template<typename Node>
inline void WithStaticValue<Type, Base>::SetValue(Node* node, Type value) {
  typedef WithValue<Nop<Type>, Base> Base1;
  assert(node);
  node->Expose()->Base1::set_value(value);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithDescValue<Group, Base>::SetValue(Node* node,
                                                 typename Group::Type value) {
  typedef WithValue<WithIdentityFilter<Group>, Base> Base1;
  Base1::SetAncDescValue(node, value);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithDescValue<Group, Base>::SetNonDescValue(
    Node* node,
    typename Group::Type value) {
  typedef WithValue<WithIdentityFilter<Group>, Base> Base1;
  node->Expose()->Base1::set_value(Group::Plus(Group::Minus(value, value),
                                               Base1::Value(node)));
}

template<typename Group, typename Base>
template<typename Node>
inline void WithAncDescValue<Group, Base>::SetValue(
    Node* node,
    typename Group::Type value) {
  typedef WithValue<Group, Base> Base1;
  Base1::SetAncDescValue(node, value);
}

template<typename Group, typename Base>
template<typename Node>
inline void WithAncDescValue<Group, Base>::SetNonDescValue(
    Node* node,
    typename Group::Type value) {
  typedef WithValue<Group, Base> Base1;
  Base1::SetValue(node,
                  Group::PlusFilter(Group::MinusFilter(value, value),
                                    Base1::Value(node)));
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type WithAggr<Semigroup, Base>::AggrAnc(
    Node* node) {
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  return AggrSeq(node);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindDirmostAnc(
    int dir,
    Node* node,
    const Predicate& predicate) {
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  return FindDirmostSeq(dir, node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindLeafmostAnc(
    Node* node,
    const Predicate& predicate) {
  return FindDirmostAnc(kLeaf, node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindRootmostAnc(
    Node* node,
    const Predicate& predicate) {
  return FindDirmostAnc(kRoot, node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type WithAggr<Semigroup, Base>::AggrProperAnc(
    Node* node) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return AggrSeq(root_cut.excl_part());
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindDirmostProperAnc(
    int dir,
    Node* node,
    const Predicate& predicate) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return FindDirmostSeq(dir, root_cut.excl_part(), predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindLeafmostProperAnc(
    Node* node,
    const Predicate& predicate) {
  return FindDirmostProperAnc(kLeaf, node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAggr<Semigroup, Base>::FindRootmostProperAnc(
    Node* node,
    const Predicate& predicate) {
  return FindDirmostProperAnc(kRoot, node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type
WithDescAggrOfDescValue<Semigroup, Base>::AggrTree(Node* node) {
  typedef WithDescAggrOfDescValue Self;
  if (!node) return Semigroup::empty_aggr();
  UpdateAggrs(node->Expose());
  return node->Self::aggr();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithDescAggrOfDescValue<Semigroup, Base>::FindDirmostTree(
    int dir,
    Node* node,
    const Predicate& predicate) {
  return (dir == kLeaf ?
          FindLeafmostTree(node, predicate) :
          FindRootmostTree(node, predicate));
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithDescAggrOfDescValue<Semigroup, Base>::FindLeafmostTree(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.FindLeafmostTree();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithDescAggrOfDescValue<Semigroup, Base>::FindRootmostTree(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.FindRootmostTree();
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type
WithDescAggrOfDescValue<Semigroup, Base>::AggrDesc(Node* node) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return AggrTree(node);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithDescAggrOfDescValue<Semigroup, Base>::FindDirmostDesc(
    int dir,
    Node* node,
    const Predicate& predicate) {
  return (dir == kLeaf ?
          FindLeafmostDesc(node, predicate) :
          FindRootmostDesc(node, predicate));
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithDescAggrOfDescValue<Semigroup, Base>::FindLeafmostDesc(
    Node* node,
    const Predicate& predicate) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return FindLeafmostTree(node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithDescAggrOfDescValue<Semigroup, Base>::FindRootmostDesc(
    Node* node,
    const Predicate& predicate) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return FindRootmostTree(node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type
WithDescAggrOfDescValue<Semigroup, Base>::AggrProperDesc(Node* node) {
  typedef WithDescAggrOfDescValue Self;
  if (!node) return Semigroup::empty_aggr();
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  ScopedCutDashed<Node> dashed_cut(node);
  if (leaf_cut.excl_part()) {
    return (dashed_cut.excl_part() ?
            Semigroup::CombineAggrs(leaf_cut.excl_part()->Self::aggr(),
                                    dashed_cut.excl_part()->Self::aggr()) :
            leaf_cut.excl_part()->Self::aggr());
  } else {
    return (dashed_cut.excl_part() ?
            dashed_cut.excl_part()->Self::aggr() :
            Semigroup::empty_aggr());
  }
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithDescAggrOfDescValue<Semigroup, Base>::FindDirmostProperDesc(
    int dir,
    Node* node,
    const Predicate& predicate) {
  return (dir == kLeaf ?
          FindLeafmostProperDesc(node, predicate) :
          FindRootmostProperDesc(node, predicate));
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithDescAggrOfDescValue<Semigroup, Base>::FindLeafmostProperDesc(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.TestProperDesc() ? search.FindLeafmostTree() : NULL;
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithDescAggrOfDescValue<Semigroup, Base>::FindRootmostProperDesc(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.TestProperDesc() ? search.FindRootmostTree() : NULL;
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
class WithDescAggrOfDescValue<Semigroup, Base>::DescSearch
    : public Search<Semigroup,
                    WithDescAggrOfDescValue<Semigroup, Base>,
                    Node,
                    Predicate> {
 public:
  DescSearch(Node* node, const Predicate& predicate)
      : Search<Semigroup, WithDescAggrOfDescValue, Node, Predicate>(
            node,
            predicate) {
  }
  Node* FindLeafmostTree() {
    while (true) {
      int i = kLeaf ^ this->node()->SolidFlipped();
      if (this->TestSolid(i)) continue;
      if (TestDashed()) continue;
      if (this->TestNode()) return this->node();
      if (this->TestSolid(i ^ 1)) continue;
      if (TestDotted(0)) continue;
      if (TestDotted(1)) continue;
      return NULL;
    }
  }
  Node* FindRootmostTree() {
    while (true) {
      int i = kRoot ^ this->node()->SolidFlipped();
      if (this->TestSolid(i)) continue;
      if (this->TestNode()) return this->node();
      if (this->TestSolid(i ^ 1)) continue;
      if (TestDashed()) continue;
      if (TestDotted(0)) continue;
      if (TestDotted(1)) continue;
      return NULL;
    }
  }
  bool TestProperDesc() {
    return (this->TestSolid(kLeaf ^ this->node()->SolidFlipped()) ||
            TestDashed());
  }
 private:
  bool TestNonSolid(Node* non_solid) {
    if (!non_solid) return false;
    non_solid->DottedAddValues(this->node());
    if (!this->TestAggr(non_solid)) {
      non_solid->DottedSubtractValues(this->node());
      return false;
    }
    this->set_node(non_solid);
    return true;
  }
  bool TestDashed() { return TestNonSolid(this->node()->dashed()); }
  bool TestDotted(int i) { return TestNonSolid(this->node()->dotted(i)); }
};

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type
WithAncDescAggr<Semigroup, Base>::AggrTree(Node* node) {
  typedef WithAncDescAggr Self;
  if (!node) return Semigroup::empty_aggr();
  UpdateAggrs(node->Expose());
  return node->Self::combined_aggr();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAncDescAggr<Semigroup, Base>::FindDirmostTree(
    int dir,
    Node* node,
    const Predicate& predicate) {
  return (dir == kLeaf ?
          FindLeafmostTree(node, predicate) :
          FindRootmostTree(node, predicate));
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithAncDescAggr<Semigroup, Base>::FindLeafmostTree(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.FindLeafmostTree();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithAncDescAggr<Semigroup, Base>::FindRootmostTree(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.FindRootmostTree();
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type
WithAncDescAggr<Semigroup, Base>::AggrDesc(Node* node) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return AggrTree(node);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAncDescAggr<Semigroup, Base>::FindDirmostDesc(
    int dir,
    Node* node,
    const Predicate& predicate) {
  return (dir == kLeaf ?
          FindLeafmostDesc(node, predicate) :
          FindRootmostDesc(node, predicate));
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAncDescAggr<Semigroup, Base>::FindLeafmostDesc(
    Node* node,
    const Predicate& predicate) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return FindLeafmostTree(node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAncDescAggr<Semigroup, Base>::FindRootmostDesc(
    Node* node,
    const Predicate& predicate) {
  ScopedCutDirOf<Node> root_cut(kRoot, node);
  return FindRootmostTree(node, predicate);
}

template<typename Semigroup, typename Base>
template<typename Node>
inline typename Semigroup::Type
WithAncDescAggr<Semigroup, Base>::AggrProperDesc(Node* node) {
  typedef WithAncDescAggr Self;
  if (!node) return Semigroup::empty_aggr();
  ScopedCutDirOf<Node> leaf_cut(kLeaf, node);
  ScopedCutDashed<Node> dashed_cut(node);
  if (leaf_cut.excl_part()) {
    return (
        dashed_cut.excl_part() ?
        Semigroup::CombineAggrs(leaf_cut.excl_part()->Self::combined_aggr(),
                                dashed_cut.excl_part()->Self::other_aggr()) :
        leaf_cut.excl_part()->Self::combined_aggr());
  } else {
    return (dashed_cut.excl_part() ?
            dashed_cut.excl_part()->Self::other_aggr() :
            Semigroup::empty_aggr());
  }
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
inline Node* WithAncDescAggr<Semigroup, Base>::FindDirmostProperDesc(
    int dir,
    Node* node,
    const Predicate& predicate) {
  return (dir == kLeaf ?
          FindLeafmostProperDesc(node, predicate) :
          FindRootmostProperDesc(node, predicate));
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithAncDescAggr<Semigroup, Base>::FindLeafmostProperDesc(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.FindLeafmostProperDesc();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
Node* WithAncDescAggr<Semigroup, Base>::FindRootmostProperDesc(
    Node* node,
    const Predicate& predicate) {
  if (!node) return NULL;
  DescSearch<Node, Predicate> search(node, predicate);
  return search.FindRootmostProperDesc();
}

template<typename Semigroup, typename Base>
template<typename Node, typename Predicate>
class WithAncDescAggr<Semigroup, Base>::DescSearch
    : public Search<Semigroup,
                    WithAncDescAggr<Semigroup, Base>,
                    Node,
                    Predicate> {
 public:
  DescSearch(Node* node, const Predicate& predicate)
      : Search<Semigroup, WithAncDescAggr, Node, Predicate>(node, predicate) {
  }
  Node* FindLeafmostTree() {
    while (true) {
      int i = kLeaf ^ this->node()->SolidFlipped();
      if (Node* solid = this->node()->solid(i)) {
        solid->SolidAddValues(this->node());
        if (DescendSolidOther(solid)) continue;
        if (TestSolid(solid)) return this->FindDirmostSeq(kLeaf);
        solid->SolidSubtractValues(this->node());
      }
      if (DescendDashed()) continue;
      if (this->TestNode()) return this->node();
      if (Node* solid = this->node()->solid(i ^ 1)) {
        solid->SolidAddValues(this->node());
        if (DescendSolidOther(solid)) continue;
        if (TestSolid(solid)) return this->FindDirmostSeq(kLeaf);
        solid->SolidSubtractValues(this->node());
      }
      return NULL;
    }
  }
  Node* FindRootmostTree() {
    while (true) {
      int i = kRoot ^ this->node()->SolidFlipped();
      if (Node* solid = this->node()->solid(i)) {
        solid->SolidAddValues(this->node());
        if (TestSolid(solid)) return this->FindDirmostSeq(kRoot);
        if (DescendSolidOther(solid)) continue;
        solid->SolidSubtractValues(this->node());
      }
      if (this->TestNode()) return this->node();
      if (Node* solid = this->node()->solid(i ^ 1)) {
        solid->SolidAddValues(this->node());
        if (TestSolid(solid)) return this->FindDirmostSeq(kRoot);
        if (DescendSolidOther(solid)) continue;
        solid->SolidSubtractValues(this->node());
      }
      if (DescendDashed()) continue;
      return NULL;
    }
  }
  Node* FindLeafmostProperDesc() {
    int i = kLeaf ^ this->node()->SolidFlipped();
    if (Node* solid = this->node()->solid(i)) {
      solid->SolidAddValues(this->node());
      if (DescendSolidOther(solid)) return FindLeafmostTree();
      if (TestSolid(solid)) return this->FindDirmostSeq(kLeaf);
      solid->SolidSubtractValues(this->node());
    }
    return DescendDashed() ? FindLeafmostTree() : NULL;
  }
  Node* FindRootmostProperDesc() {
    int i = kLeaf ^ this->node()->SolidFlipped();
    if (Node* solid = this->node()->solid(i)) {
      solid->SolidAddValues(this->node());
      if (TestSolid(solid)) return this->FindDirmostSeq(kRoot);
      if (DescendSolidOther(solid)) return FindRootmostTree();
      solid->SolidSubtractValues(this->node());
    }
    return DescendDashed() ? FindRootmostTree() : NULL;
  }
 private:
  bool TestOtherAggr(Node* child) {
    typedef WithAncDescAggr Self;
    if (!this->Test(child->Self::other_aggr())) return false;
    return true;
  }
  bool TestNonSolidOther(Node* non_solid) {
    if (!non_solid) return false;
    non_solid->DottedAddValues(this->node());
    if (!TestOtherAggr(non_solid)) {
      non_solid->DottedSubtractValues(this->node());
      return false;
    }
    this->set_node(non_solid);
    return true;
  }
  bool TestDotted(int i) { return TestNonSolidOther(this->node()->dotted(i)); }
  bool DescendDashed() {
    typedef WithAncDescAggr Self;
    if (!TestNonSolidOther(this->node()->dashed())) return false;
    while (!this->Test(this->node()->Self::filter_aggr())) {
      if (!TestDotted(0) && !TestDotted(1)) return false;
    }
    return true;
  }
  bool TestSolidOther(int i) {
    Node* solid = this->node()->solid(i);
    if (!solid) return false;
    solid->SolidAddValues(this->node());
    if (!TestOtherAggr(solid)) {
      solid->SolidSubtractValues(this->node());
      return false;
    }
    this->set_node(solid);
    return true;
  }
  bool DescendSolidOther(Node* solid) {
    if (!TestOtherAggr(solid)) return false;
    this->set_node(solid);
    while (!DescendDashed()) {
      if (!TestSolidOther(0) && !TestSolidOther(1)) return false;
    }
    return true;
  }
  bool TestSolid(Node* solid) {
    if (!this->TestAggr(solid)) return false;
    this->set_node(solid);
    return true;
  }
};

template<typename Base>
template<typename Node>
inline void WithReverseBy<Base>::Evert(Node* node) {
  if (node) ReverseSeq(node->HardExpose());
}
}  // namespace dtree
#endif  // DTREE_TREE_INL_H_
