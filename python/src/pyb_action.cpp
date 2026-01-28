/*!
        @file    pyb_action.cpp

        @brief   python binding of Action classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_action.h"
#include "bridge.h"
#include "Field/field.h"
#include "Action/action.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

void pybind_action(py::module m)
{
  //! Action classes
  //
  //  N.B. factory covers gauge actions only.
  //
  py::class_<Action>(m, "Action")
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) -> Action*
                  {
                    auto obj = Action::New(subtype);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("set_config", &Action::set_config, REDIRECT)

    .def("set_parameters", [](Action& self, py::kwargs kwargs) -> void
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

  //! Gauge Action
  //
  py::class_<Action_G_Plaq, Action>(m, "Action_G_Plaq", R"pybdoc(
Plaquette gauge action

Parameters
----------
beta : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs){
          auto obj = new Action_G_Plaq;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;
    
  py::class_<Action_G_Rectangle, Action>(m, "Action_G_Rectangle", R"pybdoc(
Rectangular gauge action

Parameters
----------
beta : float
c_plaq : float
c_rect : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs){
          auto obj = new Action_G_Rectangle;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;

    
  //! Fermion Action
  //
  py::class_<Action_F_Standard_lex, Action>(m, "Action_F_Standard_lex", R"pybdoc(
Standard fermion action using lexicographical site ordering
)pybdoc")
    .def(py::init([](py::object fopr_,
                     py::object force_,
                     py::object fprop_md_,
                     py::object fprop_h_,
                     py::kwargs kwargs) -> Action_F_Standard_lex*
                  {
                    auto obj = new PyWrap4<Action_F_Standard_lex,
                                           Fopr,
                                           Force,
                                           Fprop,
                                           Fprop>(fopr_, force_, fprop_md_, fprop_h_);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr"),
         py::arg("force"),
         py::arg("fprop_MD"),
         py::arg("fprop_H"),
         R"pybdoc(
Arguments
---------
fopr
force
fprop_MD
fprop_H
)pybdoc",
         REDIRECT
    )
    ;
  
  py::class_<Action_F_Rational, Action>(m, "Action_F_Rational", R"pybdoc(
Fermion action for rational approximation
)pybdoc")
    .def(py::init([](py::object fopr_langev_,
                     py::object fopr_h_,
                     py::object force_md_,
                     py::kwargs kwargs) -> Action_F_Rational*
                  {
                    auto obj = new PyWrap3<Action_F_Rational,
                                           Fopr,
                                           Fopr,
                                           Force>(fopr_langev_, fopr_h_, force_md_);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr_langev"),
         py::arg("fopr_H"),
         py::arg("fopr_MD"),
         R"pybdoc(
Arguments
---------
fopr_langev
fopr_H
fopr_MD
)pybdoc",
         REDIRECT
    )
    ;

  py::class_<Action_F_Ratio_lex, Action>(m, "Action_F_Ratio_lex", R"pybdoc(
Fermion action for Hasenbusch preconditioning
)pybdoc")
    .def(py::init([](py::object fopr_prec_,
                     py::object force_prec_,
                     py::object fopr_,
                     py::object force_,
                     py::object fprop_h_prec_,
                     py::object fprop_md_,
                     py::object fprop_h_,
                     py::kwargs kwargs) -> Action_F_Ratio_lex*
                  {
                    auto obj = new PyWrap7<Action_F_Ratio_lex,
                                           Fopr,
                                           Force,
                                           Fopr,
                                           Force,
                                           Fprop,
                                           Fprop,
                                           Fprop>(
                                             fopr_prec_,
                                             force_prec_,
                                             fopr_,
                                             force_,
                                             fprop_h_prec_,
                                             fprop_md_,
                                             fprop_h_);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr_prec"),
         py::arg("force_prec"),
         py::arg("fopr"),
         py::arg("force"),
         py::arg("fprop_H_prec"),
         py::arg("fprop_MD"),
         py::arg("fprop_H"),
         R"pybdoc(
Arguments
---------
fopr_prec
force_prec
fopr
force
fprop_H_prec
fprop_MD
fprop_H
)pybdoc",
         REDIRECT
    )
    ;


  //! even-odd
  //
  py::class_<Action_F_Standard_eo, Action>(m, "Action_F_Standard_eo", R"pybdoc(
Standard fermion action using even/odd site ordering
)pybdoc")
    .def(py::init([](py::object fopr_,
                     py::object force_,
                     py::object fprop_md_,
                     py::object fprop_h_,
                     py::kwargs kwargs) -> Action_F_Standard_eo*
                  {
                    auto obj = new PyWrap4<Action_F_Standard_eo,
                                           Fopr,
                                           Force,
                                           Fprop,
                                           Fprop>(fopr_, force_, fprop_md_, fprop_h_);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr"),
         py::arg("force"),
         py::arg("fprop_MD"),
         py::arg("fprop_H"),
         R"pybdoc(
Arguments
---------
fopr
force
fprop_MD
fprop_H
)pybdoc",
         REDIRECT
    )
    ;
  
  py::class_<Action_F_Ratio_eo, Action>(m, "Action_F_Ratio_eo", R"pybdoc(
Fermion action for Hasenbusch preconditioning
)pybdoc")
    .def(py::init([](py::object fopr_prec_,
                     py::object force_prec_,
                     py::object fopr_,
                     py::object force_,
                     py::object fprop_h_prec_,
                     py::object fprop_md_,
                     py::object fprop_h_,
                     py::kwargs kwargs) -> Action_F_Ratio_eo*
                  {
                    auto obj = new PyWrap7<Action_F_Ratio_eo,
                                           Fopr,
                                           Force,
                                           Fopr,
                                           Force,
                                           Fprop,
                                           Fprop,
                                           Fprop>(
                                             fopr_prec_,
                                             force_prec_,
                                             fopr_,
                                             force_,
                                             fprop_h_prec_,
                                             fprop_md_,
                                             fprop_h_);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr_prec"),
         py::arg("force_prec"),
         py::arg("fopr"),
         py::arg("force"),
         py::arg("fprop_H_prec"),
         py::arg("fprop_MD"),
         py::arg("fprop_H"),
         R"pybdoc(
Arguments
---------
fopr_prec
force_prec
fopr
force
fprop_H_prec
fprop_MD
fprop_H
)pybdoc",
         REDIRECT
    )
    ;

#ifdef TRUNK

  py::class_<Action_F_Staggered_eo, Action>(m, "Action_F_Staggered_eo", R"pybdoc(
Staggered fermion action
)pybdoc")
    .def(py::init([](py::object fopr,
                     py::object force,
                     py::kwargs kwargs) -> Action_F_Staggered_eo*
                  {
                    auto obj = new PyWrap2<Action_F_Staggered_eo,
                                           Fopr,
                                           Force>(fopr, force);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr"),
         py::arg("force"),
         R"pybdoc(
Arguments
---------
fopr
force
)pybdoc",
         REDIRECT
    )
    ;

  py::class_<Action_F_Overlap_Nf2, Action>(m, "Action_F_Overlap_Nf2", R"pybdoc(
Nf=2 Overlap fermion action
)pybdoc")
    .def(py::init([](py::object fopr,
                     py::object force,
                     py::object fprop_md,
                     py::object fprop_h,
                     py::kwargs kwargs) -> Action_F_Overlap_Nf2*
                  {
                    auto obj = new PyWrap4<Action_F_Overlap_Nf2,
                                           Fopr,
                                           Force,
                                           Fprop,
                                           Fprop>(
                                             fopr,
                                             force,
                                             fprop_md,
                                             fprop_h);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("fopr"),
         py::arg("force"),
         py::arg("fprop_MD"),
         py::arg("fprop_H"),
         R"pybdoc(
Arguments
---------
fopr
force
fprop_MD
fprop_H
)pybdoc",
         REDIRECT
    )
    ;
#endif


}
