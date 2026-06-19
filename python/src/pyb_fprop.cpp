/*!
        @file    pyb_fprop.cpp

        @brief   python binding of fermion propagator classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2021-02-27 09:53:33 #$

        @version $LastChangedRevision: 2182 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_fprop.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Fermion/fprop.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

//XXX workaround

class PyFprop_Wilson_Shift : public Fprop_Wilson_Shift
{
 public:
  PyFprop_Wilson_Shift(Fopr_Wilson *fopr)
    : Fprop_Wilson_Shift(fopr, nullptr) {}
};


void pybind_fprop(py::module m)
{
  //! Fprop
  //
  //  parameters:
  //    none
  //
  py::class_<Fprop>(m, "Fprop", R"pybdoc(
Fermion propagator class family
)pybdoc")
    .def(py::init([](const std::string& subtype, py::object solver, py::kwargs kwargs) -> Fprop*
                  {
                    Fprop *obj = nullptr;

                    if (subtype == "Standard_lex") {
                      obj = new PyWrap1<Fprop_Standard_lex, Solver>(solver);
                    }
                    
                    // if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("type"), py::arg("solver"),
         REDIRECT
    )

    .def("invert_D", [](Fprop &self, Field& x, const Field& b)
         {
           int nconv = -1;
           double diff = 0.0;
           self.invert_D(x, b, nconv, diff);
           return std::make_tuple(nconv, diff);
         }, REDIRECT)

    .def("invert_DdagD", [](Fprop &self, Field& x, const Field& b)
         {
           int nconv = -1;
           double diff = 0.0;
           self.invert_DdagD(x, b, nconv, diff);
           return std::make_tuple(nconv, diff);
         }, REDIRECT)

    .def("set_config", &Fprop::set_config, REDIRECT)
    ;

  //! Fprop_Standard_lex
  //
  //  parameters:
  //    none
  //
  py::class_<Fprop_Standard_lex, Fprop>(m, "Fprop_Standard_lex", R"pybdoc(
Calculate quark propagator for fermion operator with lexicographical site index
)pybdoc")
    .def(py::init([](py::object solver, py::kwargs kwargs){
          auto obj = new PyWrap1<Fprop_Standard_lex, Solver>(solver);
          // if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }),
      py::arg("solver"),
      REDIRECT
    )
    ;

  //! Fprop_Wilson_Shift
  //
  py::class_<Fprop_Wilson_Shift>(m, "Fprop_Wilson_Shift", R"pybdoc(
Calculate shifted quark propagator
)pybdoc")
    .def(py::init([](py::object fopr, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<PyFprop_Wilson_Shift, Fopr_Wilson>(fopr);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         REDIRECT)

    .def("invert_D", &Fprop_Wilson_Shift::invert_D, REDIRECT)
    ;

  //! Fprop_Standard_eo
  //
  //  parameters:
  //    none
  //
  py::class_<Fprop_Standard_eo, Fprop>(m, "Fprop_Standard_eo", R"pybdoc(
Calculate quark propagator for fermion operator with even-odd site index
)pybdoc")
    .def(py::init([](py::object solver, py::kwargs kwargs){
          auto obj = new PyWrap1<Fprop_Standard_eo, Solver>(solver);
          // if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }),
      py::arg("solver"),
      REDIRECT
    )
    ;

#ifdef TRUNK
  //! Fprop_Staggered_eo
  //
  //  constructor argument:
  //    Fopr_Staggered_eo fopr
  //  note:
  //    not a subclass of Fprop.
  //
  py::class_<Fprop_Staggered_eo>(m, "Fprop_Staggered_eo", R"pybdoc(
Calculate quark propagator for staggered fermion operator with even-odd site index
)pybdoc")
    .def(py::init([](py::object fopr, py::kwargs kwargs){
          auto obj = new PyWrap1<Fprop_Staggered_eo, Fopr_Staggered_eo>(fopr);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }),
      py::arg("fopr"),
      REDIRECT
    )

    .def("invert_D", &Fprop_Staggered_eo::invert_D, REDIRECT)
    ;

  //! Fprop_Domainwall_4d
  //
  //  constructor arguments:
  //    Fopr fopr_dw
  //    Solver solver
  //
  //  note:
  //    not a subclass of Fprop
  //
  py::class_<Fprop_Domainwall_4d>(m, "Fprop_Domainwall_4d", R"pybdoc(
Calculate four-dimensional space quark propagator of Domain-wall fermion
)pybdoc")
    .def(py::init([](py::object fopr_dw,
                     py::object solver,
                     py::kwargs kwargs) -> Fprop_Domainwall_4d*
                  {
                    auto obj = new PyWrap2<Fprop_Domainwall_4d, Fopr, Solver>(fopr_dw, solver);
                    // if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr"),
         py::arg("solver"),
         REDIRECT
    )

    .def("set_config", &Fprop_Domainwall_4d::set_config, REDIRECT)
    
    .def("invert_D", [](Fprop_Domainwall_4d& self, Field& dst, const Field& src)
         {
           int nconv = -1;
           double diff = 0.0;
           self.invert_D(dst, src, nconv, diff);
           return std::make_tuple(nconv, diff);
         },
         REDIRECT)
    .def("invert_DdagD", &Fprop_Domainwall_4d::invert_DdagD, REDIRECT)
    ;

  //! Fprop_Overlap_5d
  //
  //  constructor arguments:
  //    Fopr_Overlap fopr
  //    Field_G u
  //
  //  note:
  //    not a subclass of Fprop
  //
  py::class_<Fprop_Overlap_5d>(m, "Fprop_Overlap_5d", R"pybdoc(
Calculate overlap propagator with 5d solver
)pybdoc")
    .def(py::init([](py::object fopr_ov,
                     py::object u,
                     py::kwargs kwargs) -> Fprop_Overlap_5d*
                  {
                    auto obj = new PyWrap2<Fprop_Overlap_5d, Fopr_Overlap, Field_G>(fopr_ov, u);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr"),
         py::arg("u"),
         REDIRECT
    )

    .def("invert_D", &Fprop_Overlap_5d::invert_D, REDIRECT)
    
    .def("set_lowmodes", &Fprop_Overlap_5d::set_lowmodes, REDIRECT)
    ;
  
#endif
}
