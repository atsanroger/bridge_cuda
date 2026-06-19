/*!
        @file    pyb_wilsonloop.h

        @brief   python binding of Wilson loop measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_WILSONLOOP_INCLUDED
#define PYB_WILSONLOOP_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_wilsonloop(py::module m);

#endif /* PYB_WILSONLOOP_INCLUDED */
