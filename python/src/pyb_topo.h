/*!
        @file    pyb_topo.h

        @brief   python binding of topological charge measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_TOPO_INCLUDED
#define PYB_TOPO_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_topo(py::module m);

#endif /* PYB_TOPO_INCLUDED */
