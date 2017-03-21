// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef UTIL_SCOPED_H_
#define UTIL_SCOPED_H_

#include <stddef.h>

#include "util/disallow_copy_and_assign.h"

namespace util {

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
  DISALLOW_COPY_AND_ASSIGN(ScopedPtr);
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
  DISALLOW_COPY_AND_ASSIGN(ScopedArray);
};

template<typename T>
inline void swap(ScopedArray<T>& a, ScopedArray<T>& b) {
  a.swap(b);
}
}  // namespace util
#endif  // UTIL_SCOPED_H_
