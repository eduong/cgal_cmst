// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef DTREE_TYPE_H_
#define DTREE_TYPE_H_

#include <assert.h>
#include <algorithm>
#include <limits>

#define DTREE_VERSION (201405L)

namespace dtree {

template<typename Group>
struct WithTrivialFilter : public Group {
  template<typename U>
  static U PlusFilter(U x, typename Group::Type /*y*/) { return x; }
  template<typename U>
  static U MinusFilter(U x, typename Group::Type /*y*/) { return x; }
};

template<typename Group>
struct WithIdentityFilter : public Group {
  template<typename U>
  static U PlusFilter(U x, typename Group::Type y) {
    return Group::Plus(x, y);
  }
  template<typename U>
  static U MinusFilter(U x, typename Group::Type y) {
    return Group::Minus(x, y);
  }
};

template<typename Type_>
struct Nop_ {
  typedef Type_ Type;
  template<typename U>
  static U Plus(U x, Type /*y*/) { return x; }
  template<typename U>
  static U Minus(U x, Type /*y*/) { return x; }
};

template<typename Type>
struct Nop : public WithTrivialFilter<Nop_<Type> > {};

template<typename CountType, typename SumType>
struct CasAggr {
  CasAggr() : count(0), sum(0) {}
  CasAggr(CountType c, SumType s) : count(c), sum(s) {}
  CountType count;
  SumType sum;
};

template<typename Type_>
struct Add {
  typedef Type_ Type;
  template<typename U>
  static U Plus(U x, Type y) { return x + y; }
  template<typename U>
  static U Minus(U x, Type y) { return x - y; }
  template<typename CountType, typename SumType>
  static CasAggr<CountType, SumType> Plus(CasAggr<CountType, SumType> x,
                                          Type y) {
    return CasAggr<CountType, SumType>(x.count, x.sum + x.count * y);
  }
  template<typename CountType, typename SumType>
  static CasAggr<CountType, SumType> Minus(CasAggr<CountType, SumType> x,
                                           Type y) {
    return CasAggr<CountType, SumType>(x.count, x.sum - x.count * y);
  }
};

template<typename Type_>
struct Xor_ {
  typedef Type_ Type;
  template<typename U>
  static U Plus(U x, Type y) { return x ^ y; }
  template<typename U>
  static U Minus(U x, Type y) { return x ^ y; }
  static int FlippedFromValue(Type x) { return x; }
  static Type flip_delta() { return 1; }
};

template<typename Type>
struct Xor : public WithTrivialFilter<Xor_<Type> > {};

////////////////////

template<typename Type>
struct Bounds {
  static Type bottom() {
    return (std::numeric_limits<Type>::has_infinity ?
            -std::numeric_limits<Type>::infinity() :
            std::numeric_limits<Type>::min());
  }
  static Type top() {
    return (std::numeric_limits<Type>::has_infinity ?
            std::numeric_limits<Type>::infinity() :
            std::numeric_limits<Type>::max());
  }
};

template<typename Type_>
struct Min {
  typedef Type_ Type;
  template<typename U>
  static Type AggrFromValue(U x) { return x; }
  static Type CombineAggrs(Type x, Type y) { return std::min(x, y); }
  static Type empty_aggr() { return Bounds<Type>::top(); }
};

template<typename U>
class Less_ {
 public:
  explicit Less_(U bound) : bound_(bound) {}
  template<typename Type>
  bool operator()(Type x) const { return x < bound_; }
 private:
  U bound_;
};

template<typename U>
inline Less_<U> Less(U bound) {
  return Less_<U>(bound);
}

template<typename U>
class LessEqual_ {
 public:
  explicit LessEqual_(U bound) : bound_(bound) {}
  template<typename Type>
  bool operator()(Type x) const { return x <= bound_; }
 private:
  U bound_;
};

template<typename U>
inline LessEqual_<U> LessEqual(U bound) {
  return LessEqual_<U>(bound);
}

template<typename Type_>
struct Max {
  typedef Type_ Type;
  template<typename U>
  static Type AggrFromValue(U x) { return x; }
  static Type CombineAggrs(Type x, Type y) { return std::max(x, y); }
  static Type empty_aggr() { return Bounds<Type>::bottom(); }
};

template<typename U>
class Greater_ {
 public:
  explicit Greater_(U bound) : bound_(bound) {}
  template<typename Type>
  bool operator()(Type x) const { return x > bound_; }
 private:
  U bound_;
};

template<typename U>
inline Greater_<U> Greater(U bound) {
  return Greater_<U>(bound);
}

template<typename U>
class GreaterEqual_ {
 public:
  explicit GreaterEqual_(U bound) : bound_(bound) {}
  template<typename Type>
  bool operator()(Type x) const { return x >= bound_; }
 private:
  U bound_;
};

template<typename U>
inline GreaterEqual_<U> GreaterEqual(U bound) {
  return GreaterEqual_<U>(bound);
}

struct NoValue {};

template<typename Type_>
struct Count {
  typedef Type_ Type;
  static Type AggrFromValue(NoValue /*x*/) { return 1; }
  static Type CombineAggrs(Type x, Type y) { return x + y; }
  static Type empty_aggr() { return 0; }
};

template<typename Type_>
struct Sum {
  typedef Type_ Type;
  template<typename U>
  static Type AggrFromValue(U x) { return x; }
  static Type CombineAggrs(Type x, Type y) { return x + y; }
  static Type empty_aggr() { return 0; }
};

template<typename U>
class Index_ {
 public:
  explicit Index_(U i) : i_(i) { assert(i >= 0); }
  template<typename Type>
  bool operator()(Type x) const { return i_ < x; }
 private:
  U i_;
};

template<typename U>
inline Index_<U> Index(U bound) {
  return Index_<U>(bound);
}

template<typename CountType, typename SumType>
struct CountAndSum {
  typedef CasAggr<CountType, SumType> Type;
  template<typename U>
  static Type AggrFromValue(U x) { return CasAggr<CountType, SumType>(1, x); }
  static CasAggr<CountType, SumType> CombineAggrs(
      CasAggr<CountType, SumType> x,
      CasAggr<CountType, SumType> y) {
    return CasAggr<CountType, SumType>(x.count + y.count, x.sum + y.sum);
  }
  static Type empty_aggr() { return CasAggr<CountType, SumType>(); }
};

template<typename U>
class IndexByCount_ {
 public:
  explicit IndexByCount_(U i) : i_(i) { assert(i >= 0); }
  template<typename CountType, typename SumType>
  bool operator()(CasAggr<CountType, SumType> x) const { return i_ < x.count; }
 private:
  U i_;
};

template<typename U>
inline IndexByCount_<U> IndexByCount(U bound) {
  return IndexByCount_<U>(bound);
}

template<typename U>
class IndexBySum_ {
 public:
  explicit IndexBySum_(U i) : i_(i) { assert(i >= 0); }
  template<typename CountType, typename SumType>
  bool operator()(CasAggr<CountType, SumType> x) const { return i_ < x.sum; }
 private:
  U i_;
};

template<typename U>
inline IndexBySum_<U> IndexBySum(U bound) {
  return IndexBySum_<U>(bound);
}

template<typename Type_>
struct Or {
  typedef Type_ Type;
  template<typename U>
  static Type AggrFromValue(U x) { return x; }
  static Type CombineAggrs(Type x, Type y) { return x | y; }
  static Type empty_aggr() { return 0; }
};

class Nonzero_ {
 public:
  Nonzero_() {}
  template<typename Type>
  bool operator()(Type x) const { return x; }
};

inline Nonzero_ Nonzero() {
  return Nonzero_();
}

template<typename U>
class NonzeroAnd_ {
 public:
  explicit NonzeroAnd_(U mask) : mask_(mask) {}
  template<typename Type>
  bool operator()(Type x) const { return x & mask_; }
 private:
  U mask_;
};

template<typename U>
inline NonzeroAnd_<U> NonzeroAnd(U mask) {
  return NonzeroAnd_<U>(mask);
}

////////////////////////////////////////

static const int kLw = 0;
static const int kRw = 1;

template<typename Base>
class DpAggr {
 public:
  DpAggr() : dw_() {}
  DpAggr(Base lw, Base rw) : dw_() {
    set_lw(lw);
    set_rw(rw);
  }
  Base dw(int dir) const { return dw_[dir]; }
  void set_dw(int dir, Base dw_dir) { dw_[dir] = dw_dir; }
  Base lw() const { return dw(kLw); }
  void set_lw(Base lw) { set_dw(kLw, lw); }
  Base rw() const { return dw(kRw); }
  void set_rw(Base rw) { set_dw(kRw, rw); }
 private:
  Base dw_[2];
};

template<typename Base>
class DpValue : public DpAggr<Base> {
 public:
  DpValue() : flipped_(0) {}
  explicit DpValue(DpAggr<Base> x, int flipped = 0)
      : DpAggr<Base>(x), flipped_(flipped) {
  }
  DpValue(Base lw, Base rw, int flipped = 0)
      : DpAggr<Base>(lw, rw), flipped_(flipped) {
  }
  int flipped() const { return flipped_; }
  void set_flipped(int flipped) { flipped_ = flipped; }
 private:
  int flipped_;
};

template<typename Base>
class DpAdd_ {
 public:
  typedef DpValue<Base> Type;
  static DpValue<Base> Plus(DpValue<Base> x, DpValue<Base> y) {
    return DpValue<Base>(AggrPlus(x, y), x.flipped() ^ y.flipped());
  }
  static DpValue<Base> Minus(DpValue<Base> x, DpValue<Base> y) {
    return DpValue<Base>(AggrMinus(x, y), x.flipped() ^ y.flipped());
  }
  static int FlippedFromValue(DpValue<Base> x) { return x.flipped(); }
  static DpValue<Base> flip_delta() { return DpValue<Base>(0, 0, 1); }
  static DpAggr<Base> Plus(DpAggr<Base> x, DpValue<Base> y) {
    return AggrPlus(x, y);
  }
  static DpAggr<Base> Minus(DpAggr<Base> x, DpValue<Base> y) {
    return AggrMinus(x, y);
  }

 private:
  static DpAggr<Base> AggrPlus(DpAggr<Base> x, DpValue<Base> y) {
    return DpAggr<Base>(x.dw(kLw ^ y.flipped()) + y.lw(),
                        x.dw(kRw ^ y.flipped()) + y.rw());
  }
  static DpAggr<Base> AggrMinus(DpAggr<Base> x, DpValue<Base> y) {
    return DpAggr<Base>(x.dw(kLw ^ y.flipped()) - y.dw(kLw ^ y.flipped()),
                        x.dw(kRw ^ y.flipped()) - y.dw(kRw ^ y.flipped()));
  }
};

template<typename Base>
struct DpAdd : public WithTrivialFilter<DpAdd_<Base> > {};

////////////////////

template<typename Base>
struct DpMin : public DpAggr<Base> {
  typedef DpAggr<Base> Type;
  static DpAggr<Base> AggrFromValue(DpValue<Base> x) { return x; }
  static DpAggr<Base> CombineAggrs(DpAggr<Base> x, DpAggr<Base> y) {
    return DpAggr<Base>(std::min(x.lw(), y.lw()), std::min(x.rw(), y.rw()));
  }
  static DpAggr<Base> empty_aggr() {
    return DpAggr<Base>(Bounds<Base>::top(), Bounds<Base>::top());
  }
};

template<typename U>
class DwLess_ {
 public:
  DwLess_(int dir, U bound) : dir_(dir), bound_(bound) {}
  template<typename Type>
  bool operator()(DpAggr<Type> x) const { return x.dw(dir_) < bound_; }
 private:
  int dir_;
  U bound_;
};

template<typename U>
inline DwLess_<U> DwLess(int dir, U bound) {
  return DwLess_<U>(dir, bound);
}

template<typename U>
inline DwLess_<U> LwLess(U bound) {
  return DwLess_<U>(kLw, bound);
}

template<typename U>
inline DwLess_<U> RwLess(U bound) {
  return DwLess_<U>(kRw, bound);
}

template<typename U>
class DwLessEqual_ {
 public:
  DwLessEqual_(int dir, U bound) : dir_(dir), bound_(bound) {}
  template<typename Type>
  bool operator()(DpAggr<Type> x) const { return x.dw(dir_) <= bound_; }
 private:
  int dir_;
  U bound_;
};

template<typename U>
inline DwLessEqual_<U> DwLessEqual(int dir, U bound) {
  return DwLessEqual_<U>(dir, bound);
}

template<typename U>
inline DwLessEqual_<U> LwLessEqual(U bound) {
  return DwLessEqual_<U>(kLw, bound);
}

template<typename U>
inline DwLessEqual_<U> RwLessEqual(U bound) {
  return DwLessEqual_<U>(kRw, bound);
}

template<typename Base>
struct DpMax : public DpAggr<Base> {
  typedef DpAggr<Base> Type;
  static DpAggr<Base> AggrFromValue(DpValue<Base> x) { return x; }
  static DpAggr<Base> CombineAggrs(DpAggr<Base> x, DpAggr<Base> y) {
    return DpAggr<Base>(std::max(x.lw(), y.lw()), std::max(x.rw(), y.rw()));
  }
  static DpAggr<Base> empty_aggr() {
    return DpAggr<Base>(Bounds<Base>::bottom(), Bounds<Base>::bottom());
  }
};

template<typename U>
class DwGreater_ {
 public:
  DwGreater_(int dir, U bound) : dir_(dir), bound_(bound) {}
  template<typename Type>
  bool operator()(DpAggr<Type> x) const { return x.dw(dir_) > bound_; }
 private:
  int dir_;
  U bound_;
};

template<typename U>
inline DwGreater_<U> DwGreater(int dir, U bound) {
  return DwGreater_<U>(dir, bound);
}

template<typename U>
inline DwGreater_<U> LwGreater(U bound) {
  return DwGreater_<U>(kLw, bound);
}

template<typename U>
inline DwGreater_<U> RwGreater(U bound) {
  return DwGreater_<U>(kRw, bound);
}

template<typename U>
class DwGreaterEqual_ {
 public:
  DwGreaterEqual_(int dir, U bound) : dir_(dir), bound_(bound) {}
  template<typename Type>
  bool operator()(DpAggr<Type> x) const { return x.dw(dir_) >= bound_; }
 private:
  int dir_;
  U bound_;
};

template<typename U>
inline DwGreaterEqual_<U> DwGreaterEqual(int dir, U bound) {
  return DwGreaterEqual_<U>(dir, bound);
}

template<typename U>
inline DwGreaterEqual_<U> LwGreaterEqual(U bound) {
  return DwGreaterEqual_<U>(kLw, bound);
}

template<typename U>
inline DwGreaterEqual_<U> RwGreaterEqual(U bound) {
  return DwGreaterEqual_<U>(kRw, bound);
}
}  // namespace dtree
#endif  // DTREE_TYPE_H_
