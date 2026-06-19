/*!
        @file    pyb_gauge_fixing.cpp

        @brief   python binding of gauge fixing classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_gauge_fixing.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Gauge/gaugeFixing.h"
#include "Parameters/parameters.h"

void pybind_gauge_fixing(py::module m)
{
  //! GaugeFixing
  //
  py::class_<GaugeFixing>(m, "GaugeFixing", R"pybdoc(
Performs gauge fixing

See Also
--------
GaugeFixing_Landau
GaugeFixing_Coulomb
)pybdoc")
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) -> GaugeFixing*
                  {
                    auto obj = GaugeFixing::New(subtype);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("fix", &GaugeFixing::fix, REDIRECT)

    .def("set_parameters", [](GaugeFixing& self, py::kwargs kwargs){
        self.set_parameters(make_parameters(kwargs));
      }, REDIRECT)
    ;

  py::class_<GaugeFixing_Landau, GaugeFixing>(m, "GaugeFixing_Landau", R"pybdoc(
Landau gauge fixing

Parameters
----------
maximum_number_of_iteration : int
number_of_naive_iteration : int
interval_of_measurement : int
iteration_to_reset : int
convergence_criterion_squared : float
overrelaxation_parameter : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> GaugeFixing_Landau*
                  {
                    auto obj = new GaugeFixing_Landau;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<GaugeFixing_Coulomb, GaugeFixing>(m, "GaugeFixing_Coulomb", R"pybdoc(
Coulomb gauge fixing

Parameters
----------
maximum_number_of_iteration : int
number_of_naive_iteration : int
interval_of_measurement : int
iteration_to_reset : int
convergence_criterion_squared : float
overrelaxation_parameter : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> GaugeFixing_Coulomb*
                  {
                    auto obj = new GaugeFixing_Coulomb;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

}
