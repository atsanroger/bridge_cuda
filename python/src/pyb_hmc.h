/*!
        @file    pyb_hmc.h

        @brief   python binding of HMC classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_HMC_INCLUDED
#define PYB_HMC_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_hmc(py::module m);

#endif /* PYB_HMC_INCLUDED */
