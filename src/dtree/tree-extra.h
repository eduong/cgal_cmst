// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef DTREE_TREE_EXTRA_H_
#define DTREE_TREE_EXTRA_H_


namespace dtree {

template<typename T>
class ScopedPtr {
 public:
  explicit ScopedPtr(T* ptr) : ptr_(ptr) {}
  ~ScopedPtr() { delete ptr_; }
  void reset(T* ptr) { ScopedPtr(ptr).swap(*this); }
  T& operator*() const { return *ptr_; }
  T* operator->() const { return ptr_; }
  T* get() const { return ptr_; }
  void swap(ScopedPtr& b) {
    T* tmp = ptr_;
    ptr_ = b.ptr_;
    b.ptr_ = tmp;
  }
 private:
  T* ptr_;
  ScopedPtr(const ScopedPtr&);
  void operator=(const ScopedPtr&);
};

template<typename T>
inline void swap(ScopedPtr<T>& a, ScopedPtr<T>& b) {
  a.swap(b);
}

template<typename T>
class ScopedArray {
 public:
  explicit ScopedArray(T* array) : array_(array) {}
  ~ScopedArray() { delete[] array_; }
  void reset(T* array) { ScopedArray(array).swap(*this); }
  T& operator[](ptrdiff_t i) const { return array_[i]; }
  T* get() const { return array_; }
  void swap(ScopedArray& b) {
    T* tmp = array_;
    array_ = b.array_;
    b.array_ = tmp;
  }
 private:
  T* array_;
  ScopedArray(const ScopedArray&);
  void operator=(const ScopedArray&);
};

template<typename T>
inline void swap(ScopedArray<T>& a, ScopedArray<T>& b) {
  a.swap(b);
}

template<typename Node, typename Visitor>
void EvertByTraversing(Node* node, Visitor* visitor) {
  if (!node) return;
  Node* u = node->HardExpose();
  while (true) {
    while (Node* c = u->solid(kLeaf ^ u->SolidFlipped())) {
      c->SolidAddValues(u);
      u = c;
    }
    Node* e = u->solid(kRoot ^ u->SolidFlipped());
    if (e) {
      e->SolidAddValues(u);
      (*visitor)(u);
      u->set_solid(kLeaf ^ u->SolidFlipped(), NULL);
      u->set_solid(kRoot ^ u->SolidFlipped(), e);
      u = e;
      continue;
    }
    (*visitor)(u);
    u->SolidResetAggrs();
    while (true) {
      Node* p = u->parent();
      if (!p) return;
      Node* w = p->solid(kRoot ^ p->SolidFlipped());
      if (w == u) {
        w = p->solid(kLeaf ^ p->SolidFlipped());
        p->SolidResetAggrs();
        p->set_solid(kLeaf ^ p->SolidFlipped(), u);
        u->SolidSubtractValues(p);
        p->SolidMergeSolidAggrs(u);
        p->set_solid(kRoot ^ p->SolidFlipped(), w);
        if (w) {
          w->SolidSubtractValues(p);
          p->SolidMergeSolidAggrs(w);
        }
        u = p;
        continue;
      }
      if (w) {
        w->SolidAddValues(p);
        (*visitor)(p);
        p->set_solid(kLeaf ^ p->SolidFlipped(), u);
        p->set_solid(kRoot ^ p->SolidFlipped(), w);
        u = w;
        break;
      }
      (*visitor)(p);
      p->SolidResetAggrs();
      p->set_solid(kRoot ^ p->SolidFlipped(), u);
      u->SolidSubtractValues(p);
      p->SolidMergeSolidAggrs(u);
      p->set_solid(kLeaf ^ p->SolidFlipped(), NULL);
      u = p;
    }
  }
}

template<typename Node,
         typename Prim,
         void (Node::*Visit)(Node* path[]),
         void (Node::*VisitNonRoot)()>
Node* AssembleHelper(Node* path[], size_t path_size, size_t zero) {
  if (path_size == 1) {
    Node* u = path[0];
    (u->*Visit)(path + path_size);
    Prim::ResetAggrs(u);
    Prim::set_child(u, 0, NULL);
    Prim::set_child(u, 1, NULL);
    return u;
  }
  size_t key = zero + (path[0]->num_desc() - zero) / 2;
  size_t i = 0;
  size_t n = path_size;
  while (true) {
    if (n == 1) break;
    size_t h = n / 2;
    size_t j = i + h;
    if (path[j]->num_desc() > key) {
      i = j;
      n -= h;
    } else {
      n = h;
    }
  }
  Node* u = path[i];
  size_t u_num_desc = u->num_desc();
  (u->*Visit)(path + path_size);
  Prim::ResetAggrs(u);
  size_t j = i + 1;
  Node* c;
  if (j < path_size) {
    c = AssembleHelper<Node, Prim, Visit, VisitNonRoot>(path + j,
                                                        path_size - j, zero);
    (c->*VisitNonRoot)();
    Prim::SubtractValues(c, u);
    Prim::MergeAggrs(u, c);
    c->set_parent(u);
  } else {
    c = NULL;
  }
  Prim::set_child(u, kLeaf ^ Prim::Flipped(u), c);
  Node* e;
  if (i > 0) {
    e = AssembleHelper<Node, Prim, Visit, VisitNonRoot>(path, i, u_num_desc);
    (e->*VisitNonRoot)();
    Prim::SubtractValues(e, u);
    Prim::MergeAggrs(u, e);
    e->set_parent(u);
  } else {
    e = NULL;
  }
  Prim::set_child(u, kRoot ^ Prim::Flipped(u), e);
  return u;
}

template<typename Base>
inline void EndTree<Base>::ResetDashed() {
}

template<typename Base>
inline void EndTree<Base>::ResetDotted() {
}

template<typename Base>
inline void EndTree<Base>::AssembleDashed(EndTree* /*path*/[]) {
}

template<typename Base>
inline void EndTree<Base>::AssembleSolidRoot(EndTree* parent) {
  set_parent(parent);
}

template<typename Base>
inline void EndTreeWithDesc<Base>::ResetDashed() {
  set_dashed(NULL);
}

template<typename Base>
inline void EndTreeWithDesc<Base>::ResetDotted() {
  set_dotted(0, NULL);
  set_dotted(1, NULL);
}

template<typename Base>
inline void EndTreeWithDesc<Base>::AssembleSolidRoot(
    EndTreeWithDesc* parent) {
  if (parent) {
    set_heavy(parent->dashed());
    parent->set_dashed(this);
  } else {
    set_parent(NULL);
    ResetDotted();
  }
}

template<typename Base>
inline void EndTreeWithDesc<Base>::AssembleDashed(EndTreeWithDesc* path[]) {
  typedef EndTreeWithDesc Node;
  Node* u = dashed();
  if (!u) return;
  size_t path_size = 0;
  size_t cum_num_desc = 0;
  do {
    path[path_size++] = u;
    cum_num_desc += u->num_desc();
    u->set_num_desc(cum_num_desc);
    u = u->heavy();
  } while (u);
  std::reverse(path, path + path_size);
  typedef DottedPrimitives Prim;
  u = AssembleHelper<Node, Prim, &Node::SolidToDashedAggrs, &Node::nop>(
      path, path_size, 0);
  set_dashed(u);
  u->set_parent(this);
}

template<typename Base>
inline void EndTreeWithDesc<Base>::SolidToDashedAggrs(
    EndTreeWithDesc* /*path*/[]) {
  Base::SolidToDashedAggrs();
}

template<typename Node>
inline void AssembleTopologicallySorted(Node* scratch[], size_t num_nodes) {
  Node** solid_root = scratch;
  for (Node** i = scratch; i != scratch + num_nodes; ++i) {
    Node* u = *i;
    Node* p = u->parent();
    if (!p) {
      *solid_root++ = u;
      continue;
    }
    p->set_num_desc(p->num_desc() + u->num_desc());
    Node* h = p->heavy();
    if (!h) {
      p->set_heavy(u);
    } else if (u->num_desc() > h->num_desc()) {
      *solid_root++ = h;
      p->set_heavy(u);
    } else {
      *solid_root++ = u;
    }
  }
  std::reverse(scratch, solid_root);
  while (solid_root != scratch) {
    Node** path = --solid_root;
    size_t path_size = 1;
    Node* u = *path;
    Node* p = u->parent();
    while (true) {
      u = u->heavy();
      if (!u) break;
      path[path_size++] = u;
    }
    typedef SolidPrimitives<Node> Prim;
    u = AssembleHelper<Node, Prim, &Node::AssembleDashed, &Node::ResetDotted>(
        path, path_size, 0);
    if (p) u->DottedSubtractValues(p);
    u->AssembleSolidRoot(p);
  }
}

template<typename IteratorLike>
void Assemble(IteratorLike first, IteratorLike last) {
  size_t num_nodes = 0;
  typedef typename IteratorLike::value_type NodePointer;
  for (IteratorLike i = first; i != last; ++i) {
    NodePointer u = *i;
    u->ResetDashed();
    u->set_num_desc(0);
    u->set_heavy(NULL);
    ++num_nodes;
  }
  ScopedArray<NodePointer> scratch(new NodePointer[num_nodes]);
  NodePointer* j = scratch.get() + num_nodes;
  for (IteratorLike i = first; i != last; ++i) {
    NodePointer u = *i;
    if (u->num_desc() != 0) continue;
    NodePointer* k = j;
    do {
      u->set_num_desc(1);
      *--j = u;
      u = u->parent();
    } while (u && u->num_desc() == 0);
    std::reverse(j, k);
  }
  AssembleTopologicallySorted(scratch.get(), num_nodes);
}

template<typename T>
class Contiguous {
 public:
  typedef T value_type;
  Contiguous() : value_() {}
  explicit Contiguous(T value) : value_(value) {}
  bool operator!=(Contiguous x) const { return value_ != x.value_; }
  T operator*() const { return value_; }
  Contiguous& operator++() {
    ++value_;
    return *this;
  }
 private:
  T value_;
};

template<typename Node>
inline void Assemble(Node node[], size_t num_nodes) {
  Assemble(Contiguous<Node*>(node), Contiguous<Node*>(node + num_nodes));
}

template<typename Node>
void CutOneOfMany(Node* node) {
  Node* u = node;
  u->ResetDashed();
  u->ResetDotted();
  Node* p = u->parent();
  if (!p) return;
  Node* c = NULL;
  do {
    u->set_parent(c);
    c = u;
    u = p;
    p = u->parent();
  } while (p);
  while (true) {
    int i;
    if (c->SolidIndex(u, &i)) {
      u->set_solid(i, NULL);
      c->SolidAddValues(u);
    } else {
      c->DottedAddValues(u);
    }
    u = c;
    c = u->parent();
    if (!c) break;
    u->set_parent(NULL);
  }
}
}  // namespace dtree
#endif  // DTREE_TREE_EXTRA_H_
