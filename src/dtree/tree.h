// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef DTREE_TREE_H_
#define DTREE_TREE_H_

#include "seq.h"

namespace dtree {

template<typename Base>
class EndTree : public Base {
 public:
  EndTree() : parent_(NULL), solid_() {}
  EndTree* parent() const { return parent_; }
  void set_parent(EndTree* parent) { parent_ = parent; }
  EndTree* solid(int i) const { return solid_[i]; }
  void set_solid(int i, EndTree* solid_i) { solid_[i] = solid_i; }
  template<typename Forwarder>
  void CopyFrom(const EndTree& node, Forwarder forward);
  void SolidResetAggrs() { Base::SolidResetAggrs(NULL); }
  bool SolidIndex(const EndTree* parent, int* i) const {
    if (!parent) return false;
    if (parent->solid(0) == this) {
      *i = 0;
      return true;
    }
    if (parent->solid(1) == this) {
      *i = 1;
      return true;
    }
    return false;
  }
  template<bool finishing_search>
  void ExposeStep(EndTree** d0, EndTree** d1);
  template<bool finishing_search>
  void DashedReplace(EndTree* c, EndTree* d0, EndTree* d1);
  template<bool finishing_search>
  void DashedRemove(EndTree* d0, EndTree* d1);
  void DashedInsert(EndTree* c);
  EndTree* Expose();
  EndTree* HardExpose();
  void FinishSearch();

 public:
  size_t num_desc() const { return abstract_.num_desc; }
  void set_num_desc(size_t num_desc) { abstract_.num_desc = num_desc; }
  EndTree* heavy() const { return abstract_.heavy; }
  void set_heavy(EndTree* heavy) { abstract_.heavy = heavy; }
  void ResetDashed();
  void ResetDotted();
  void AssembleDashed(EndTree* path[]);
  void AssembleSolidRoot(EndTree* parent);

 private:
  EndTree* parent_;
  union {
    EndTree* solid_[2];
    struct {
      size_t num_desc;
      EndTree* heavy;
    } abstract_;
  };
  EndTree(const EndTree&);
  void operator=(const EndTree&);
};

template<typename Base>
class EndTreeWithDesc : public Base {
 public:
  struct DottedPrimitives;
  EndTreeWithDesc() : parent_(NULL), dashed_(NULL), solid_(), dotted_() {}
  EndTreeWithDesc* parent() const { return parent_; }
  void set_parent(EndTreeWithDesc* parent) { parent_ = parent; }
  EndTreeWithDesc* dashed() const { return dashed_; }
  void set_dashed(EndTreeWithDesc* dashed) { dashed_ = dashed; }
  EndTreeWithDesc* solid(int i) const { return solid_[i]; }
  void set_solid(int i, EndTreeWithDesc* solid_i) { solid_[i] = solid_i; }
  EndTreeWithDesc* dotted(int i) const { return dotted_[i]; }
  void set_dotted(int i, EndTreeWithDesc* dotted_i) { dotted_[i] = dotted_i; }
  template<typename Forwarder>
  void CopyFrom(const EndTreeWithDesc& node, Forwarder forward);
  void SolidResetAggrs() { Base::SolidResetAggrs(dashed()); }
  bool SolidIndex(const EndTreeWithDesc* parent, int* i) const {
    if (!parent) return false;
    if (parent->solid(0) == this) {
      *i = 0;
      return true;
    }
    if (parent->solid(1) == this) {
      *i = 1;
      return true;
    }
    return false;
  }
  void SolidToDashed(const EndTreeWithDesc* parent) {
    Base::SolidToDashed(parent);
    Base::SolidToDashedAggrs();
  }
  void DottedResetAggrs() {
    Base::DottedResetAggrs(dashed());
    if (EndTreeWithDesc* c = solid(0)) this->DottedMergeSolidAggrs(c);
    if (EndTreeWithDesc* e = solid(1)) this->DottedMergeSolidAggrs(e);
  }
  bool DottedIndex(const EndTreeWithDesc* parent, int* i) const {
    if (parent->dotted(0) == this) {
      *i = 0;
      return true;
    }
    if (parent->dotted(1) == this) {
      *i = 1;
      return true;
    }
    return false;
  }
  template<bool finishing_search>
  void ExposeStep(EndTreeWithDesc** d0, EndTreeWithDesc** d1);
  template<bool finishing_search>
  void DashedReplace(EndTreeWithDesc* c,
                     EndTreeWithDesc* d0,
                     EndTreeWithDesc* d1);
  template<bool finishing_search>
  void DashedRemove(EndTreeWithDesc* d0, EndTreeWithDesc* d1);
  void DashedInsert(EndTreeWithDesc* c);
  EndTreeWithDesc* Expose();
  EndTreeWithDesc* HardExpose();
  void FinishSearch();
  EndTreeWithDesc* DottedCut(int i);
  void DottedLink(int i, EndTreeWithDesc* dotted_i);
  template<bool finishing_search>
  void DottedSplay(EndTreeWithDesc** d0, EndTreeWithDesc** d1);
  void Splice(EndTreeWithDesc* d);
  void AbsorbDashed(EndTreeWithDesc* v);

 public:
  size_t num_desc() const { return abstract_.num_desc; }
  void set_num_desc(size_t num_desc) { abstract_.num_desc = num_desc; }
  EndTreeWithDesc* heavy() const { return abstract_.heavy; }
  void set_heavy(EndTreeWithDesc* heavy) { abstract_.heavy = heavy; }
  void ResetDashed();
  void ResetDotted();
  void AssembleDashed(EndTreeWithDesc* path[]);
  void AssembleSolidRoot(EndTreeWithDesc* parent);
  void SolidToDashedAggrs(EndTreeWithDesc* path[]);
  void nop() {}

 private:
  EndTreeWithDesc* parent_;
  EndTreeWithDesc* dashed_;
  EndTreeWithDesc* solid_[2];
  union {
    EndTreeWithDesc* dotted_[2];
    struct {
      size_t num_desc;
      EndTreeWithDesc* heavy;
    } abstract_;
  };
  EndTreeWithDesc(const EndTreeWithDesc&);
  void operator=(const EndTreeWithDesc&);
};

template<typename Type, typename Base>
class WithStaticValue : public WithValue<Nop<Type>, Base> {
 public:
  enum { kAncUpdates = false };
  template<typename Node>
  static Type Value(Node* node);
  template<typename Node>
  static void SetValue(Node* node, Type value);
 protected:
  WithStaticValue() {}
};

#define DTREE_BASE1 WithValue<WithTrivialFilter<Group>, Base>
template<typename Group, typename Base>
class WithAncValue : public DTREE_BASE1 {
 public:
  using DTREE_BASE1::AddToAnc;
  using DTREE_BASE1::SubtractFromAnc;
  using DTREE_BASE1::AddToProperAnc;
  using DTREE_BASE1::SubtractFromProperAnc;
 protected:
  WithAncValue() {}
};
#undef DTREE_BASE1

#define DTREE_BASE1 WithValue<WithIdentityFilter<Group>, Base>
template<typename Group, typename Base>
class WithDescValue : public DTREE_BASE1 {
 public:
  enum { kAncUpdates = false };
  template<typename Node>
  static void SetValue(Node* node, typename Group::Type value);
  template<typename Node>
  static void SetNonDescValue(Node* node, typename Group::Type value);
  using DTREE_BASE1::AddToTree;
  using DTREE_BASE1::SubtractFromTree;
  using DTREE_BASE1::AddToDesc;
  using DTREE_BASE1::SubtractFromDesc;
  using DTREE_BASE1::AddToProperDesc;
  using DTREE_BASE1::SubtractFromProperDesc;
 protected:
  WithDescValue() {}
};
#undef DTREE_BASE1

#define DTREE_BASE1 WithValue<Group, Base>
template<typename Group, typename Base>
class WithAncDescValue : public DTREE_BASE1 {
 public:
  template<typename Node>
  static void SetValue(Node* node, typename Group::Type value);
  template<typename Node>
  static void SetNonDescValue(Node* node, typename Group::Type value);
  using DTREE_BASE1::AddToAnc;
  using DTREE_BASE1::SubtractFromAnc;
  using DTREE_BASE1::AddToProperAnc;
  using DTREE_BASE1::SubtractFromProperAnc;
  using DTREE_BASE1::AddToTree;
  using DTREE_BASE1::SubtractFromTree;
  using DTREE_BASE1::AddToDesc;
  using DTREE_BASE1::SubtractFromDesc;
  using DTREE_BASE1::AddToProperDesc;
  using DTREE_BASE1::SubtractFromProperDesc;
 protected:
  WithAncDescValue() {}
};
#undef DTREE_BASE1

#define DTREE_BASE1 WithAggr<Semigroup, Base>
template<typename Semigroup, typename Base>
class WithAncAggr : public DTREE_BASE1 {
 public:
  using DTREE_BASE1::AggrAnc;
  using DTREE_BASE1::FindDirmostAnc;
  using DTREE_BASE1::FindLeafmostAnc;
  using DTREE_BASE1::FindRootmostAnc;
  using DTREE_BASE1::AggrProperAnc;
  using DTREE_BASE1::FindDirmostProperAnc;
  using DTREE_BASE1::FindLeafmostProperAnc;
  using DTREE_BASE1::FindRootmostProperAnc;
 protected:
  WithAncAggr() {}
};
#undef DTREE_BASE1

template<typename Semigroup, typename Base>
class WithDescAggrOfDescValue : public Base {
 private:
  typedef typename Base::Group_ Group;

 public:
  typedef Semigroup Semigroup_;
  typename Semigroup::Type delta_aggr() const { return delta_aggr_; }
  void set_delta_aggr(typename Semigroup::Type delta_aggr) {
    delta_aggr_ = delta_aggr;
  }
  typename Semigroup::Type singleton_aggr() const {
    return Semigroup::AggrFromValue(this->value());
  }
  typename Semigroup::Type aggr() const {
    return Group::Plus(delta_aggr(), this->value());
  }
  void SolidResetAggrs(const WithDescAggrOfDescValue* dashed) {
    Base::SolidResetAggrs(dashed);
    Reset1(dashed);
  }
  void SolidMergeSolidAggrs(const WithDescAggrOfDescValue* solid) {
    Base::SolidMergeSolidAggrs(solid);
    Merge1(solid);
  }
  void DottedResetAggrs(const WithDescAggrOfDescValue* dashed) {
    Base::DottedResetAggrs(dashed);
    Reset1(dashed);
  }
  void DottedMergeSolidAggrs(const WithDescAggrOfDescValue* solid) {
    Base::DottedMergeSolidAggrs(solid);
    Merge1(solid);
  }
  void DottedMergeDottedAggrs(const WithDescAggrOfDescValue* dotted) {
    Base::DottedMergeDottedAggrs(dotted);
    Merge1(dotted);
  }

 public:
  template<typename Node>
  static typename Semigroup::Type AggrTree(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostTree(int dir,
                               Node* node,
                               const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostTree(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostTree(Node* node, const Predicate& predicate);
  template<typename Node>
  static typename Semigroup::Type AggrDesc(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostDesc(int dir,
                               Node* node,
                               const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostDesc(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostDesc(Node* node, const Predicate& predicate);
  template<typename Node>
  static typename Semigroup::Type AggrProperDesc(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostProperDesc(int dir,
                                     Node* node,
                                     const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostProperDesc(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostProperDesc(Node* node, const Predicate& predicate);

 protected:
  WithDescAggrOfDescValue() : delta_aggr_() {}

 private:
  template<typename Node, typename Predicate>
  class DescSearch;
  void Reset1(const WithDescAggrOfDescValue* dashed) {
    set_delta_aggr(Group::Minus(singleton_aggr(), this->value()));
    if (dashed) Merge1(dashed);
  }
  void Merge1(const WithDescAggrOfDescValue* child) {
    set_delta_aggr(Semigroup::CombineAggrs(delta_aggr(), child->aggr()));
  }
  typename Semigroup::Type delta_aggr_;

 private:
  void AggrSeq();
  void FindDirmostSeq();
  void FindLeftmostSeq();
  void FindRightmostSeq();
  void AggrAnc();
  void FindDirmostAnc();
  void FindLeafmostAnc();
  void FindRootmostAnc();
  void AggrProperAnc();
  void FindDirmostProperAnc();
  void FindLeafmostProperAnc();
  void FindRootmostProperAnc();
};

#define DTREE_BASE1 WithAggr<Semigroup, Base>
template<typename Semigroup, typename Base>
class WithAncDescAggr : public DTREE_BASE1 {
 private:
  typedef typename DTREE_BASE1::Group_ Group;

 public:
  typename Semigroup::Type delta_filter_aggr() const {
    return Group::Plus(this->delta_aggr(),
                       Group::MinusFilter(this->value(), this->value()));
  }
  typename Semigroup::Type filter_aggr() const {
    return Group::PlusFilter(this->delta_aggr(), this->value());
  }
  typename Semigroup::Type delta_other_aggr() const {
    return delta_other_aggr_;
  }
  void set_delta_other_aggr(typename Semigroup::Type delta_other_aggr) {
    delta_other_aggr_ = delta_other_aggr;
  }
  typename Semigroup::Type other_aggr() const {
    return Group::PlusFilter(delta_other_aggr(), this->value());
  }
  typename Semigroup::Type combined_aggr() const {
    return Semigroup::CombineAggrs(this->aggr(), other_aggr());
  }
  void SolidResetAggrs(const WithAncDescAggr* dashed) {
    DTREE_BASE1::SolidResetAggrs(dashed);
    set_delta_other_aggr(dashed ?
                         dashed->other_aggr() :
                         Semigroup::empty_aggr());
  }
  void SolidMergeSolidAggrs(const WithAncDescAggr* solid) {
    DTREE_BASE1::SolidMergeSolidAggrs(solid);
    Merge1(solid);
  }
  void SolidToDashedAggrs() {
    DTREE_BASE1::SolidToDashedAggrs();
    this->set_delta_aggr(Semigroup::CombineAggrs(delta_filter_aggr(),
                                                 delta_other_aggr()));
  }
  void DottedResetAggrs(const WithAncDescAggr* dashed) {
    DTREE_BASE1::DottedResetAggrs(dashed);
    set_delta_other_aggr(this->delta_aggr());
  }
  void DottedMergeDottedAggrs(const WithAncDescAggr* dotted) {
    DTREE_BASE1::DottedMergeDottedAggrs(dotted);
    Merge1(dotted);
  }

 public:
  using DTREE_BASE1::AggrAnc;
  using DTREE_BASE1::FindDirmostAnc;
  using DTREE_BASE1::FindLeafmostAnc;
  using DTREE_BASE1::FindRootmostAnc;
  using DTREE_BASE1::AggrProperAnc;
  using DTREE_BASE1::FindDirmostProperAnc;
  using DTREE_BASE1::FindLeafmostProperAnc;
  using DTREE_BASE1::FindRootmostProperAnc;
  template<typename Node>
  static typename Semigroup::Type AggrTree(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostTree(int dir,
                               Node* node,
                               const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostTree(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostTree(Node* node, const Predicate& predicate);
  template<typename Node>
  static typename Semigroup::Type AggrDesc(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostDesc(int dir,
                               Node* node,
                               const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostDesc(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostDesc(Node* node, const Predicate& predicate);
  template<typename Node>
  static typename Semigroup::Type AggrProperDesc(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostProperDesc(int dir,
                                     Node* node,
                                     const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostProperDesc(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostProperDesc(Node* node, const Predicate& predicate);

 protected:
  WithAncDescAggr() : delta_other_aggr_(Semigroup::empty_aggr()) {}

 private:
  template<typename Node, typename Predicate>
  class DescSearch;
  void Merge1(const WithAncDescAggr* child) {
    set_delta_other_aggr(Semigroup::CombineAggrs(delta_other_aggr(),
                                                 child->other_aggr()));
  }
  typename Semigroup::Type delta_other_aggr_;
};
#undef DTREE_BASE1

template<bool anc_updates>
struct WdaSelector;

template<>
struct WdaSelector<false> {
  template<typename Semigroup, typename Base>
  class Wda : public WithDescAggrOfDescValue<Semigroup, Base> {
   protected:
    Wda() {}
  };
};

template<>
struct WdaSelector<true> {
  template<typename Semigroup, typename Base>
  class Wda : public WithAncDescAggr<Semigroup, Base> {
   protected:
    Wda() {}
  };
};

template<typename Semigroup, typename Base>
class WithDescAggr
    : public WdaSelector<Base::kAncUpdates>::template Wda<Semigroup, Base> {
 protected:
  WithDescAggr() {}
};

#define DTREE_BASE1 WithReverseBy<Base>
template<typename Base>
class WithEvertBy : public DTREE_BASE1 {
 public:
  using DTREE_BASE1::Evert;
 protected:
  WithEvertBy() {}
};
#undef DTREE_BASE1

template<typename Base>
class WithEvert : public WithEvertBy<WithAncValue<Xor<int>, Base> > {
 protected:
  WithEvert() {}
};
}  // namespace dtree
#endif  // DTREE_TREE_H_
