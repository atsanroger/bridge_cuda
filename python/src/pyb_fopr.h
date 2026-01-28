/*!
        @file    pyb_fopr.h

        @brief   python binding of Fermion operator classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_FOPR_INCLUDED
#define PYB_FOPR_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_fopr(py::module m);

#endif /* PYB_FOPR_INCLUDED */
