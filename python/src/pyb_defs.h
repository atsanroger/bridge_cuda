/*!
        @file    pyb_defs.h

        @brief   definitions of utility functions

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_DEFS_INCLUDED
#define PYB_DEFS_INCLUDED

#include <iostream>
#include <cstdio>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/iostream.h>
namespace py = pybind11;

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef DEBUG
#define ENTER do { std::cout << "+++ " << __PRETTY_FUNCTION__ << " [" << this << "]" << std::endl; } while(0)
#define ENTER_F do { std::cout << "+++ " << __PRETTY_FUNCTION__ << std::endl; } while(0)
#else
#define ENTER
#define ENTER_F
#endif

#define warn(fmt, ...)                                  \
  do {                                                  \
    char buf[1024];                                     \
    snprintf(buf, 1024, "Warn: " fmt, ##__VA_ARGS__);   \
    std::cout << buf << std::endl;                      \
  } while(0)


class __attribute__((visibility("hidden")))
custom_redirector : public py::scoped_ostream_redirect
{
 private:
  using base = py::scoped_ostream_redirect;
 public:
  custom_redirector()
    : base(vout.getStream(), py::module::import("sys").attr("stdout"))
  {}
};

#define REDIRECT py::call_guard<custom_redirector>()

#endif /* PYB_DEFS_INCLUDED */
