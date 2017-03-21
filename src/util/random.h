// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef UTIL_RANDOM_H_
#define UTIL_RANDOM_H_

#include <assert.h>
#include <stdint.h>

#include "util/disallow_copy_and_assign.h"

namespace util {

static const uint64_t kDefaultSeed = UINT64_C(88172645463325252);

// Marsaglia, "Xorshift RNGs", xor64()
class Random {
 public:
  Random() : seed_(kDefaultSeed) {}
  uint64_t seed() const { return seed_; }
  void set_seed(uint64_t seed) { seed_ = seed ? seed : kDefaultSeed; }
  uint64_t Next() {
    seed_ ^= seed() << 13;
    seed_ ^= seed() >> 7;
    seed_ ^= seed() << 17;
    return seed();
  }
  template<typename T>
  T Next(T n) {
    assert(n > 0);
    uint64_t upper_bound = UINT64_MAX - (UINT64_MAX - (n - 1)) % n;
    while (Next() > upper_bound) {}
    return seed() % n;
  }
 private:
  uint64_t seed_;
  DISALLOW_COPY_AND_ASSIGN(Random);
};
}  // namespace util
#endif  // UTIL_RANDOM_H_
