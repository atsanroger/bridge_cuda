/*!
        @file    pyb_shift_field.h

        @brief   python binding of field shift classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_SHIFT_FIELD_INCLUDED
#define PYB_SHIFT_FIELD_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_shift_field(py::module m);

#endif /* PYB_SHIFT_FIELD_INCLUDED */
