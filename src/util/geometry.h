// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef UTIL_GEOMETRY_H_
#define UTIL_GEOMETRY_H_

#include <limits>

namespace util {

struct GridPoint {
  GridPoint()
      : x(std::numeric_limits<int>::min()),
        y(std::numeric_limits<int>::min()) {
  }
  GridPoint(int x_, int y_) : x(x_), y(y_) {}
  int x;
  int y;
};

struct BoundingBox {
  BoundingBox()
      : xmin(std::numeric_limits<double>::quiet_NaN()),
        ymin(std::numeric_limits<double>::quiet_NaN()),
        xmax(std::numeric_limits<double>::quiet_NaN()),
        ymax(std::numeric_limits<double>::quiet_NaN()) {
  }
  BoundingBox(double xmin_, double ymin_, double xmax_, double ymax_)
      : xmin(xmin_), ymin(ymin_), xmax(xmax_), ymax(ymax_) {
  }
  double xmin;
  double ymin;
  double xmax;
  double ymax;
};

struct Point {
  Point()
      : x(std::numeric_limits<double>::quiet_NaN()),
        y(std::numeric_limits<double>::quiet_NaN()) {
  }
  Point(double x_, double y_) : x(x_), y(y_) {}
  double x;
  double y;
};
}  // namespace util
#endif  // UTIL_GEOMETRY_H_
