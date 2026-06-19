/*!
        @file    pyb_rng.h

        @brief   python binding of random number generator classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_RNG_INCLUDED
#define PYB_RNG_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_rng(py::module m);

#endif /* PYB_IO_INCLUDED */
