/*!
        @file    pyb_wilsonloop.cpp

        @brief   python binding of Wilson loop measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_wilsonloop.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Gauge/wilsonLoop.h"
#include "Measurements/Gauge/polyakovLoop.h"
#include "Parameters/parameters.h"

void pybind_wilsonloop(py::module m)
{
  //! WilsonLoop
  //
  //  parameters:
  //    max_spatial_loop_size (int)
  //    max_temporal_loop_size (int)
  //    number_of_loop_type (int)
  //
  py::class_<WilsonLoop>(m, "WilsonLoop", R"pybdoc(
Wilson loop measurement

Parameters
----------
max_spatial_loop_size : int
max_temporal_loop_size : int
number_of_loop_type : int
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs){
          auto obj = new WilsonLoop;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("measure", &WilsonLoop::measure, REDIRECT)
    ;

  //! PolyakovLoop
  //
  //  parameters:
  //    spatial_correlator_size (int)
  //    number_of_correlator_type (int)
  //
  py::class_<PolyakovLoop>(m, "PolyakovLoop", R"pybdoc(
Polyakov loop measurement

Paramters
---------
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs){
          auto obj = new PolyakovLoop;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("measure", &PolyakovLoop::measure_ploop, REDIRECT)
    // .def("measure_correlator", &PolyakovLoop::measure_ploop_corr, REDIRECT)
    ;
          
}
