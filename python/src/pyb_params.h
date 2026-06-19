/*!
        @file    pyb_params.h

        @brief   translating python dictionary to Parameters class

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_PARAMS_INCLUDED
#define PYB_PARAMS_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "Parameters/parameters.h"

Parameters make_parameters(py::dict kwargs);

void pybind_params(py::module& m);

#endif /* PYB_PARAMS_INCLUDED */
