/*!
        @file    pyb_common.h

        @brief   interface to CommonParameters

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_COMMON_INCLUDED
#define PYB_COMMON_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_common(py::module& m);

#endif /* PYB_COMMON_INCLUDED */
