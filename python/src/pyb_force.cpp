/*!
        @file    pyb_force.cpp

        @brief   python binding of force calculation classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_force.h"
#include "bridge.h"
#include "Field/field.h"
#include "Force/Gauge/force_G.h"
#include "Force/Fermion/force_F.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

void pybind_force(py::module m)
{
  //! Gauge Force
  //
  py::class_<Force_G>(m, "Force_G", R"pybdoc(
Gauge force calculation
)pybdoc")
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) -> Force_G*
                  {
                    auto obj = Force_G::New(subtype);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("set_config", (void (Force_G::*)(Field_G *)) &Force_G::set_config, REDIRECT)
    ;

  py::class_<Force_G_Plaq, Force_G>(m, "Force_G_Plaq", R"pybdoc(
Gauge force calculation for plaquette gauge action

Parameters
----------
beta : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs){
          auto obj = new Force_G_Plaq;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;

  py::class_<Force_G_Rectangle, Force_G>(m, "Force_G_Rectangle", R"pybdoc(
Gauge force calculation for rectangular gauge action

Parameters
----------
beta : float
c_plaq : float
c_rect : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs){
          auto obj = new Force_G_Rectangle;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;

  //! Fermion Force
  //
  py::class_<Force>(m, "Force", R"pybdoc(
Fermion force calculation

See Also
--------
Force_F_Wilson_Nf2
Force_F_Clover_Nf2
Force_F_Wilson_Nf2_Isochemical
Force_F_Clover_Nf2_Isochemical
Force_F_Wilson_eo
)pybdoc"
#ifdef TRUNK
R"pybdoc(
Force_F_Wilson_TwistedMass_Nf2
Force_F_Staggered_eo
Force_F_Domainwall
Force_F_Overlap_Nf2
)pybdoc"
#endif
R"pybdoc(
Force_F_Smeared
Force_F_Rational
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Force*
                  {
                    warn("invalid instantiation of fermion force base class");
                    return nullptr;
                  }), REDIRECT)

    .def("set_config", (void (Force::*)(Field *)) &Force::set_config, REDIRECT)
    ;

  py::class_<Force_F_Wilson_Nf2, Force>(m, "Force_F_Wilson_Nf2")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Wilson_Nf2*
                  {
                    std::string repr = kwargs.contains("gamma_matrix_type")
                      ? kwargs["gamma_matrix_type"].cast<std::string>()
                      : "Dirac";

                    auto obj = new Force_F_Wilson_Nf2(repr);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
  py::class_<Force_F_Clover_Nf2, Force>(m, "Force_F_Clover_Nf2")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Clover_Nf2*
                  {
                    std::string repr = kwargs.contains("gamma_matrix_type")
                      ? kwargs["gamma_matrix_type"].cast<std::string>()
                      : "Dirac";

                    auto obj = new Force_F_Clover_Nf2(repr);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
  py::class_<Force_F_Wilson_Nf2_Isochemical, Force>(m, "Force_F_Wilson_Nf2_Isochemical")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Wilson_Nf2_Isochemical*
                  {
                    std::string repr = kwargs.contains("gamma_matrix_type")
                      ? kwargs["gamma_matrix_type"].cast<std::string>()
                      : "Dirac";

                    auto obj = new Force_F_Wilson_Nf2_Isochemical(repr);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
  py::class_<Force_F_Clover_Nf2_Isochemical, Force>(m, "Force_F_Clover_Nf2_Isochemical")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Clover_Nf2_Isochemical*
                  {
                    std::string repr = kwargs.contains("gamma_matrix_type")
                      ? kwargs["gamma_matrix_type"].cast<std::string>()
                      : "Dirac";

                    auto obj = new Force_F_Clover_Nf2_Isochemical(repr);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
  py::class_<Force_F_Smeared, Force>(m, "Force_F_Smeared")
    .def(py::init([](py::object force,
                     py::object force_smear,
                     py::object director_smear,
                     py::kwargs kwargs) -> Force_F_Smeared*
                  {
                    auto obj = new PyWrap3<Force_F_Smeared,
                                           Force,
                                           ForceSmear,
                                           Director_Smear>(
                                             force,
                                             force_smear,
                                             director_smear);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Force_F_Rational, Force>(m, "Force_F_Rational")
    .def(py::init([](py::object fopr,
                     py::object force,
                     py::kwargs kwargs) -> Force_F_Rational*
                  {
                    auto obj = new PyWrap2<Force_F_Rational,
                                           Fopr,
                                           Force>(
                                             fopr,
                                             force);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;


  //! Force of smeared operator
  //
  py::class_<ForceSmear>(m, "ForceSmear")
    .def(py::init([](const std::string& subtype,
                     py::object proj, py::kwargs kwargs)
                  {
                    ForceSmear *obj = nullptr;

                    if (subtype == "APE") {
                      obj = new PyWrap1<ForceSmear_APE, Projection>(proj);
                    } else if (subtype == "HYP") {
                      obj = new PyWrap1<ForceSmear_HYP, Projection>(proj);
                    } else {
                      ;
                    }

                    if (obj) obj->set_parameters(make_parameters(kwargs));

                    return obj;
                  }), REDIRECT)
    ;

  py::class_<ForceSmear_APE, ForceSmear>(m, "ForceSmear_APE")
    .def(py::init([](py::object proj, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<ForceSmear_APE, Projection>(proj);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
  py::class_<ForceSmear_HYP, ForceSmear>(m, "ForceSmear_HYP")
    .def(py::init([](py::object proj, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<ForceSmear_HYP, Projection>(proj);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  //! Force of even-odd preconditioned fermion operators
  //
  py::class_<Force_F_Wilson_eo, Force>(m, "Force_F_Wilson_eo")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Wilson_eo*
                  {
                    std::string repr = kwargs.contains("gamma_matrix_type")
                      ? kwargs["gamma_matrix_type"].cast<std::string>()
                      : "Dirac";

                    auto obj = new Force_F_Wilson_eo(repr);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
#ifdef TRUNK

  py::class_<Force_F_Staggered_eo, Force>(m, "Force_F_Staggered_eo")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Staggered_eo*
                  {
                    auto obj = new Force_F_Staggered_eo;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Force_F_Wilson_TwistedMass_Nf2, Force>(m, "Force_F_Wilson_TwistedMass_Nf2")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Wilson_TwistedMass_Nf2*
                  {
                    std::string repr = kwargs.contains("gamma_matrix_type")
                      ? kwargs["gamma_matrix_type"].cast<std::string>()
                      : "Dirac";

                    auto obj = new Force_F_Wilson_TwistedMass_Nf2(repr);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Force_F_Domainwall, Force>(m, "Force_F_Domainwall")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Domainwall*
                  {
                    auto obj = new Force_F_Domainwall;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
  py::class_<Force_F_Overlap_Nf2, Force>(m, "Force_F_Overlap_Nf2")
    .def(py::init([](py::kwargs kwargs) -> Force_F_Overlap_Nf2*
                  {
                    auto obj = new Force_F_Overlap_Nf2;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;
  
#endif

}
