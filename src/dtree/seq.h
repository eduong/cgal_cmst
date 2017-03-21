// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef DTREE_SEQ_H_
#define DTREE_SEQ_H_

#include <stddef.h>

#include "type.h"

namespace dtree {

static const int kDesc = 0;
static const int kAnc = 1;

class Empty {
 protected:
  Empty() {}
};

template<typename Base = Empty>
class Begin : public Base {
 public:
  typedef Nop<NoValue> Group_;
  NoValue value() const { return NoValue(); }
  void SolidAddValues(const Begin* /*parent*/) {}
  void SolidSubtractValues(const Begin* /*parent*/) {}
  void SolidResetAggrs(const Begin* /*dashed*/) {}
  void SolidMergeSolidAggrs(const Begin* /*solid*/) {}
  int SolidFlipped() const { return 0; }
 public:
  enum { kAncUpdates = false };
  void SolidToDashed(const Begin* /*parent*/) {}
  void SolidToDashedAggrs() {}
  void DottedAddValues(const Begin* /*parent*/) {}
  void DottedSubtractValues(const Begin* /*parent*/) {}
  void DottedResetAggrs(const Begin* /*dashed*/) {}
  void DottedMergeSolidAggrs(const Begin* /*solid*/) {}
  void DottedMergeDottedAggrs(const Begin* /*dotted*/) {}
 protected:
  Begin() {}
};

template<typename Base>
class EndSeq : public Base {
 public:
  EndSeq() : parent_(NULL), solid_() {}
  EndSeq* parent() const { return parent_; }
  void set_parent(EndSeq* parent) { parent_ = parent; }
  EndSeq* solid(int i) const { return solid_[i]; }
  void set_solid(int i, EndSeq* solid_i) { solid_[i] = solid_i; }
  template<typename Forwarder>
  void CopyFrom(const EndSeq& node, Forwarder forward);
  void SolidResetAggrs() { Base::SolidResetAggrs(NULL); }
  bool SolidIndex(const EndSeq* parent, int* i) const {
    if (!parent) return false;
    *i = parent->solid(0) != this;
    return true;
  }
  EndSeq* Expose();
  void FinishSearch();
 private:
  EndSeq* parent_;
  EndSeq* solid_[2];
  EndSeq(const EndSeq&);
  void operator=(const EndSeq&);
};

template<typename Group, typename Base>
class WithValue : public Base {
 public:
  typedef Group Group_;
  typename Group::Type value() const { return value_; }
  void set_value(typename Group::Type value) { value_ = value; }
  void add_to_value(typename Group::Type delta) {
    set_value(Group::Plus(value(), delta));
  }
  void subtract_from_value(typename Group::Type delta) {
    set_value(Group::Minus(value(), delta));
  }
  void add_to_value_filter(typename Group::Type delta) {
    set_value(Group::PlusFilter(value(), delta));
  }
  void subtract_from_value_filter(typename Group::Type delta) {
    set_value(Group::MinusFilter(value(), delta));
  }
  void SolidAddValues(const WithValue* parent) {
    Base::SolidAddValues(parent);
    add_to_value(parent->value());
  }
  void SolidSubtractValues(const WithValue* parent) {
    Base::SolidSubtractValues(parent);
    subtract_from_value(parent->value());
  }
  void SolidToDashed(const WithValue* parent) {
    Base::SolidToDashed(parent);
    add_to_value(Group::MinusFilter(parent->value(), parent->value()));
  }
  void DottedAddValues(const WithValue* parent) {
    Base::DottedAddValues(parent);
    add_to_value_filter(parent->value());
  }
  void DottedSubtractValues(const WithValue* parent) {
    Base::DottedSubtractValues(parent);
    subtract_from_value_filter(parent->value());
  }

 public:
  template<typename Node>
  static typename Group::Type Value(Node* node);
  template<typename Node>
  static void SetValue(Node* node, typename Group::Type value);
  template<typename Node>
  static void AddToSeq(Node* node, typename Group::Type delta);
  template<typename Node>
  static void SubtractFromSeq(Node* node, typename Group::Type delta);

 public:
  enum { kAncUpdates = true };

 protected:
  template<typename Node>
  static void AddToAnc(Node* node, typename Group::Type delta);
  template<typename Node>
  static void SubtractFromAnc(Node* node, typename Group::Type delta);
  template<typename Node>
  static void AddToProperAnc(Node* node, typename Group::Type delta);
  template<typename Node>
  static void SubtractFromProperAnc(Node* node, typename Group::Type delta);

 protected:
  template<typename Node>
  static void SetAncDescValue(Node* node, typename Group::Type value);
  template<typename Node>
  static void AddToTree(Node* node, typename Group::Type delta);
  template<typename Node>
  static void SubtractFromTree(Node* node, typename Group::Type delta);
  template<typename Node>
  static void AddToDesc(Node* node, typename Group::Type delta);
  template<typename Node>
  static void SubtractFromDesc(Node* node, typename Group::Type delta);
  template<typename Node>
  static void AddToProperDesc(Node* node, typename Group::Type delta);
  template<typename Node>
  static void SubtractFromProperDesc(Node* node, typename Group::Type delta);

 protected:
  WithValue() : value_() {}

 private:
  typename Group::Type value_;

 private:
  void SetNonDescValue();
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
  void AggrTree();
  void FindDirmostTree();
  void FindLeafmostTree();
  void FindRootmostTree();
  void AggrDesc();
  void FindDirmostDesc();
  void FindLeafmostDesc();
  void FindRootmostDesc();
  void AggrProperDesc();
  void FindDirmostProperDesc();
  void FindLeafmostProperDesc();
  void FindRootmostProperDesc();
  void ReverseSeq();
  void Evert();
};

template<typename Semigroup, typename Base>
class WithAggr : public Base {
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
  void SolidResetAggrs(const WithAggr* dashed) {
    Base::SolidResetAggrs(dashed);
    set_delta_aggr(Group::Minus(singleton_aggr(), this->value()));
  }
  void SolidMergeSolidAggrs(const WithAggr* solid) {
    Base::SolidMergeSolidAggrs(solid);
    set_delta_aggr(Semigroup::CombineAggrs(delta_aggr(), solid->aggr()));
  }

 public:
  template<typename Node>
  static typename Semigroup::Type AggrSeq(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostSeq(int dir, Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeftmostSeq(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRightmostSeq(Node* node, const Predicate& predicate);

 protected:
  template<typename Node>
  static typename Semigroup::Type AggrAnc(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostAnc(int dir, Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostAnc(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostAnc(Node* node, const Predicate& predicate);
  template<typename Node>
  static typename Semigroup::Type AggrProperAnc(Node* node);
  template<typename Node, typename Predicate>
  static Node* FindDirmostProperAnc(int dir,
                                    Node* node,
                                    const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindLeafmostProperAnc(Node* node, const Predicate& predicate);
  template<typename Node, typename Predicate>
  static Node* FindRootmostProperAnc(Node* node, const Predicate& predicate);

 protected:
  WithAggr() : delta_aggr_() {}

 private:
  typename Semigroup::Type delta_aggr_;

 private:
  void AggrTree();
  void FindDirmostTree();
  void FindLeafmostTree();
  void FindRootmostTree();
  void AggrDesc();
  void FindDirmostDesc();
  void FindLeafmostDesc();
  void FindRootmostDesc();
  void AggrProperDesc();
  void FindDirmostProperDesc();
  void FindLeafmostProperDesc();
  void FindRootmostProperDesc();
};

template<typename Base>
class WithReverseBy : public Base {
 private:
  typedef typename Base::Group_ Group;
 public:
  int SolidFlipped() const {
    return Base::SolidFlipped() ^ Group::FlippedFromValue(this->value());
  }
 public:
  template<typename Node>
  static void ReverseSeq(Node* node);
 protected:
  template<typename Node>
  static void Evert(Node* node);
 protected:
  WithReverseBy() {}
};

template<typename Base>
class WithReverse : public WithReverseBy<WithValue<Xor<int>, Base> > {
 protected:
  WithReverse() {}
};
}  // namespace dtree
#endif  // DTREE_SEQ_H_
