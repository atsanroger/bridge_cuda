/*!
        @file    pyb_gradientflow.h

        @brief   python binding of gradient flow classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_GRADIENTFLOW_INCLUDED
#define PYB_GRADIENTFLOW_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_gradientflow(py::module m);

#endif /* PYB_GRADIENTFLOW_INCLUDED */
