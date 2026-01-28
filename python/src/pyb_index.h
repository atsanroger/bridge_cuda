/*!
        @file    pyb_index.h

        @brief   python binding of site index classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_INDEX_INCLUDED
#define PYB_INDEX_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_index(py::module m);

#endif /* PYB_INDEX_INCLUDED */
