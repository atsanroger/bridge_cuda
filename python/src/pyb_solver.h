/*!
        @file    pyb_solver.h

        @brief   python binding of linear equation solver classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_SOLVER_INCLUDED
#define PYB_SOLVER_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_solver(py::module m);

#endif /* PYB_SOLVER_INCLUDED */
