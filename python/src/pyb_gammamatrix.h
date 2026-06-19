/*!
        @file    pyb_gammamatrix.h

        @brief   python binding of gamma matrices

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_GAMMAMATRIX_INCLUDED
#define PYB_GAMMAMATRIX_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_gammamatrix(py::module m);

#endif /* PYB_GAMMAMATRIX_INCLUDED */
