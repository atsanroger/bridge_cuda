/*!
        @file    pyb_tools.h

        @brief   python binding of tools and mathematical functions

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_TOOLS_INCLUDED
#define PYB_TOOLS_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_tools(py::module m);

void pybind_module_tools(py::module m);

#endif /* PYB_TOOLS_INCLUDED */
