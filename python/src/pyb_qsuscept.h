/*!
        @file    pyb_qsuscept.h

        @brief   python binding of quark number susceptibility measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_QSUSCEPT_INCLUDED
#define PYB_QSUSCEPT_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_measure_qsuscept(py::module m);

#endif /* PYB_QSUSCEPT_INCLUDED */
