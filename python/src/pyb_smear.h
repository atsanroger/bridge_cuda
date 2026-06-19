/*!
        @file    pyb_smear.h

        @brief   python binding of smearing classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_SMEAR_INCLUDED
#define PYB_SMEAR_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_smear(py::module m);

#endif /* PYB_SMEAR_INCLUDED */
