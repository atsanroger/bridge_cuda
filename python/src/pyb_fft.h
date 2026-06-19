/*!
        @file    pyb_fft.h

        @brief   python binding of FFT classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYB_FFT_INCLUDED
#define PYB_FFT_INCLUDED

#include <pybind11/pybind11.h>
namespace py = pybind11;

void pybind_fft(py::module m);

#endif /* PYB_FFT_INCLUDED */
