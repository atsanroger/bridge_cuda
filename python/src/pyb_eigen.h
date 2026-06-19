/*!
        @file    pyb_eigen.h

        @brief   python binding of Eigensolver classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_EIGEN_INCLUDED
#define PYB_EIGEN_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_eigen(py::module m);

#endif /* PYB_EIGEN_INCLUDED */
