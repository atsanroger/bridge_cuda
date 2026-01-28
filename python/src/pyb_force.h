/*!
        @file    pyb_force.h

        @brief   python binding of force calculation classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_FORCE_INCLUDED
#define PYB_FORCE_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_force(py::module m);

#endif /* PYB_FORCE_INCLUDED */
