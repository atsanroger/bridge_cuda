/*!
        @file    pyb_shift_field.cpp

        @brief   python binding of field shift classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_field.h"
#include "pyb_shift_field.h"
#include "bridge.h"
#include "Field/field.h"
#include "Field/shiftField_lex.h"
#include "Parameters/parameters.h"

void pybind_shift_field(py::module m)
{
  //! ShiftField
  //
  py::class_<ShiftField_lex>(m, "ShiftField_lex", R"pybdoc(
Shift a field in lexicographical site index
)pybdoc")
    .def(py::init<>(), REDIRECT)

    .def("forward", [](ShiftField_lex& self,
                       Field& dst, const Field& src, int mu, int bc) -> void
         {
           self.forward(dst, src, bc, mu);  // note order of args
           return;
         },
         py::arg("dst"), py::arg("src"), py::arg("mu"), py::arg("boundary_condition") = 1,
         REDIRECT
    )
    
    .def("backward", [](ShiftField_lex& self,
                        Field& dst, const Field& src, int mu, int bc) -> void
         {
           self.backward(dst, src, bc, mu);  // note order of args
           return;
         },
         py::arg("dst"), py::arg("src"), py::arg("mu"), py::arg("boundary_condition") = 1,
         REDIRECT
    )
    ;
          
  //! ShiftField_eo
  //

  // enum Direction {
  //   odd_to_even = 0,
  //   even_to_odd = 1
  // };

  // py::enum_<Direction>(m, "Direction")
  //   .value("odd_to_even", Direction::odd_to_even)
  //   .value("even_to_odd", Direction::even_to_odd)
  //   .export_values();
  
  py::class_<ShiftField_eo>(m, "ShiftField_eo", R"pybdoc(
Shift a field in even-odd site index
)pybdoc")
    .def(py::init<>(), REDIRECT)

    // direction, ieo (0: e<-o, 1: o<-e)
    .def("forward_h", [](ShiftField_eo& self,
                         Field& dst, const Field& src, int mu, int ieo, int bc) -> void
         {
           self.forward_h(dst, src, bc, mu, ieo);  // note order of args
           return;
         },
         py::arg("dst"), py::arg("src"), py::arg("mu"), py::arg("ieo"), py::arg("boundary_condition") = 1,
         REDIRECT)
    
    .def("backward_h", [](ShiftField_eo& self,
                          Field& dst, const Field& src, int mu, int ieo, int bc) -> void
         {
           self.backward_h(dst, src, bc, mu, ieo);  // note order of args
           return;
         },
         py::arg("dst"), py::arg("src"), py::arg("mu"), py::arg("ieo"), py::arg("boundary_condition") = 1,
         REDIRECT)
    
    .def("forward", [](ShiftField_eo& self,
                       Field& dst, const Field& src, int mu, int bc) -> void
         {
           self.forward(dst, src, bc, mu);  // note order of args
           return;
         },
         py::arg("dst"), py::arg("src"), py::arg("mu"), py::arg("boundary_condition") = 1,
         REDIRECT
    )
    
    .def("backward", [](ShiftField_eo& self,
                        Field& dst, const Field& src, int mu, int bc) -> void
         {
           self.backward(dst, src, bc, mu);  // note order of args
           return;
         },
         py::arg("dst"), py::arg("src"), py::arg("mu"), py::arg("boundary_condition") = 1,
         REDIRECT
    )
    ;
}
