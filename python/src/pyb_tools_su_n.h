/*!
        @file    pyb_tools_su_n.h

        @brief   python binding of SU(N) matrix and vector classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_TOOLS_SU_N_INCLUDED
#define PYB_TOOLS_SU_N_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_tools_su_n(py::module m);

#endif /* PYB_TOOLS_SU_N_INCLUDED */
