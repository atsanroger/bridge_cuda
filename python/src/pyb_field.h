/*!
        @file    pyb_field.h

        @brief   python binding of Field classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_FIELD_INCLUDED
#define PYB_FIELD_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_field(py::module m);

#endif /* PYB_FIELD_INCLUDED */
