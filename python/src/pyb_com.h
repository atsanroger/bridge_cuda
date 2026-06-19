/*!
        @file    pyb_com.h

        @brief   interface to Communicator

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_COM_INCLUDED
#define PYB_COM_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_com(py::module m);

#endif /* PYB_COM_INCLUDED */
