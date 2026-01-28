/*!
        @file    pyb_staple.h

        @brief   python binding of staple and plaquette measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_STAPLE_INCLUDED
#define PYB_STAPLE_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_staple(py::module m);

#endif /* PYB_STAPLE_INCLUDED */
