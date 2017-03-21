// -*-C++-*-
// Copyright (c) 2012--2014 David Eisenstat <eisenstatdavid@gmail.com>
//     and Brown University
// Released under http://opensource.org/licenses/MIT
// May 2014 version

#ifndef UTIL_DISALLOW_COPY_AND_ASSIGN_H_
#define UTIL_DISALLOW_COPY_AND_ASSIGN_H_

#ifndef DISALLOW_COPY_AND_ASSIGN
// Google, "C++ Style Guide"
#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)
#endif
#endif  // UTIL_DISALLOW_COPY_AND_ASSIGN_H_
