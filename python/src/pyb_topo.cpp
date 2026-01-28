/*!
        @file    pyb_topo.cpp

        @brief   python binding of topological charge measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_topo.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Gauge/topologicalCharge.h"
#include "Measurements/Gauge/fieldStrength.h"
#include "Parameters/parameters.h"

void pybind_topo(py::module m)
{
  //! Topological Charge Density
  //
  py::class_<TopologicalCharge>(m, "TopologicalCharge", R"pybdoc(
Topological charge measurement

Parameters
----------
c_plaq : float
c_rect : float
max_momentum : int
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new TopologicalCharge;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

#define DEFINE_METHOD(x) def(#x, &TopologicalCharge:: x, REDIRECT)

    .DEFINE_METHOD(measure)
    .DEFINE_METHOD(measure_topological_charge)
    .DEFINE_METHOD(measure_topological_density_t)
    .DEFINE_METHOD(measure_topological_density_t_FT)
    .DEFINE_METHOD(measure_topological_density_x)
    .DEFINE_METHOD(measure_topological_density_x_FT)
    .DEFINE_METHOD(measure_topological_density_y)
    .DEFINE_METHOD(measure_topological_density_y_FT)
    .DEFINE_METHOD(measure_topological_density_z)
    .DEFINE_METHOD(measure_topological_density_z_FT)

    .DEFINE_METHOD(set_field_strength)

#undef DEFINE_METHOD

    .def("set_parameters", [](TopologicalCharge& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

  //! Field Strength
  //
  py::class_<FieldStrength>(m, "FieldStrength", R"pybdoc(
Field strength construction
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new FieldStrength;
          // if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("1x1", &FieldStrength::construct_Fmunu_1x1, REDIRECT)
    .def("1x2", &FieldStrength::construct_Fmunu_1x2, REDIRECT)
    .def("1x1_traceless", &FieldStrength::construct_Fmunu_1x1_traceless, REDIRECT)
    .def("1x2_traceless", &FieldStrength::construct_Fmunu_1x2_traceless, REDIRECT)
    ;

  //! Energy Density
  //
  py::class_<EnergyDensity>(m, "EnergyDensity", R"pybdoc(
Energy density measurement

Parameters
----------
c_plaq : float
c_rect : float
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new EnergyDensity;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("E_plaq", &EnergyDensity::E_plaq, REDIRECT)
    .def("E_clover", &EnergyDensity::E_clover, REDIRECT)

    .def("set_parameters", [](EnergyDensity& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;
    
  //! Energy Momentum Tensor
  //
  py::class_<EnergyMomentumTensor>(m, "EnergyMomentumTensor", R"pybdoc(
Energy momentum tensor measurement

Parameters
----------
c_plaq : float
c_rect : float
max_momentum : int
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new EnergyMomentumTensor;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

#define DEFINE_METHOD(x) def(#x, &EnergyMomentumTensor:: x, REDIRECT)

    .DEFINE_METHOD(measure_EMT)
    .DEFINE_METHOD(measure_EMT_at_t)
    .DEFINE_METHOD(measure_EMT_at_t_FT)
    .DEFINE_METHOD(measure_EMT_at_x)
    .DEFINE_METHOD(measure_EMT_at_x_FT)
    .DEFINE_METHOD(measure_EMT_at_y)
    .DEFINE_METHOD(measure_EMT_at_y_FT)
    .DEFINE_METHOD(measure_EMT_at_z)
    .DEFINE_METHOD(measure_EMT_at_z_FT)

    .DEFINE_METHOD(set_field_strength)
    
#undef DEFINE_METHOD

    .def("set_parameters", [](EnergyMomentumTensor& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

}
