/*!
        @file    pyb_staple.cpp

        @brief   python binding of staple and plaquette measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_staple.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Gauge/staple.h"
#include "Parameters/parameters.h"

void pybind_staple(py::module m)
{
  //! Staple and Plaquette
  //
  py::class_<Staple>(m, "Staple", R"pybdoc(
Staple construction class

See Also
--------
Staple_lex
Staple_eo
)pybdoc")
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) {
          auto obj = Staple::New(subtype);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("staple", &Staple::staple, REDIRECT)
    .def("staple_upper", &Staple::upper, REDIRECT)
    .def("staple_lower", &Staple::lower, REDIRECT)

    .def("plaquette", &Staple::plaquette, REDIRECT)
    ;

  py::class_<Staple_lex, Staple>(m, "Staple_lex", R"pybdoc(
Staple construction in lexicographical site index

Parameters
----------
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new Staple_lex;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;
  
  py::class_<Staple_eo, Staple>(m, "Staple_eo", R"pybdoc(
Staple construction in even-odd site index

Parameters
----------
filename_output : string, optional
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new Staple_eo;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;
  
}
