/*!
        @file    pyb_source.h

        @brief   python binding of source classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_SOURCE_INCLUDED
#define PYB_SOURCE_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_source(py::module m);

#endif /* PYB_SOURCE_INCLUDED */
