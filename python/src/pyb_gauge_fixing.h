/*!
        @file    pyb_gauge_fixing.h

        @brief   python binding of gauge fixing classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_GAUGE_FIXING_INCLUDED
#define PYB_GAUGE_FIXING_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_gauge_fixing(py::module m);

#endif /* PYB_GAUGE_FIXING_INCLUDED */
