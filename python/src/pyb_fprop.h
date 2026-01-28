/*!
        @file    pyb_fprop.h

        @brief   python binding of fermion propagator classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2021-02-27 09:53:33 #$

        @version $LastChangedRevision: 2182 $
*/

#ifndef PYB_FPROP_INCLUDED
#define PYB_FPROP_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_fprop(py::module m);

#endif /* PYB_FPROP_INCLUDED */
