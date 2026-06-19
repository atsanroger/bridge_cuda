/*!
        @file    pyb_io.h

        @brief   python binding of file I/O classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_IO_INCLUDED
#define PYB_IO_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_io(py::module m);

#endif /* PYB_IO_INCLUDED */
