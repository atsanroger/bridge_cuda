/*!
        @file    pyb_smear.cpp

        @brief   python binding of smearing classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <memory>
#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_smear.h"
#include "bridge.h"
#include "Smear/smear.h"
#include "Smear/projection.h"
#include "Smear/director_Smear.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

void pybind_smear(py::module m)
{
  //! Smear class
  //
  //   constructor
  //     Smear(Projection *proj)
  //
  py::class_<Smear>(m, "Smear", R"pybdoc(
Smearing of link variables

See Also
--------
Smear_APE
Smear_APE_spatial
Smear_HYP
)pybdoc")
    .def(py::init([](const std::string& subtype, py::object proj, py::kwargs kwargs) -> Smear*
                  {
                    Smear *obj = nullptr;
                    
                    if (subtype == "APE") {
                      obj = new PyWrap1<Smear_APE, Projection>(proj);
                    } else if (subtype == "APE_spatial") {
                      obj = new PyWrap1<Smear_APE_spatial, Projection>(proj);
                    } else if (subtype == "HYP") {
                      obj = new PyWrap1<Smear_HYP, Projection>(proj);
                    } else {
                      ;
                    }

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;

                  }), REDIRECT)

    .def("smear", &Smear::smear, REDIRECT)

    .def("set_parameters", [](Smear& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;

  py::class_<Smear_APE, Smear>(m, "Smear_APE", R"pybdoc(
APE-type smearing of link variables

Parameters
----------
rho_uniform : float
)pybdoc")
    .def(py::init([](py::object proj, py::kwargs kwargs) -> Smear_APE*
                  {
                    auto obj = new PyWrap1<Smear_APE, Projection>(proj);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Smear_APE_spatial, Smear>(m, "Smear_APE_spatial", R"pybdoc(
APE-type smaring of spatial link variables

Parameters
----------
rho : float
)pybdoc")
    .def(py::init([](py::object proj, py::kwargs kwargs) -> Smear_APE_spatial*
                  {
                    auto obj = new PyWrap1<Smear_APE_spatial, Projection>(proj);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Smear_HYP, Smear>(m, "Smear_HYP", R"pybdoc(
HYP smearing of link variables

Parameters
----------
alpha1 : float
alpha2 : float
alpha3 : float
)pybdoc")
    .def(py::init([](py::object proj, py::kwargs kwargs) -> Smear_HYP*
                  {
                    auto obj = new PyWrap1<Smear_HYP, Projection>(proj);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  //! Projection class
  //
  //    nominal binding. internal use.
  //
  py::class_<Projection>(m, "Projection", R"pybdoc(
Projection operator into gauge group

See Also
--------
Projection_Maximum_SU_N
Projection_Stout_SU3
)pybdoc")
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) -> Projection*
                  {
                    auto obj = Projection::New(subtype);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Projection_Maximum_SU_N, Projection>(m, "Projection_Maximum_SU_N", R"pybdoc(
Maximum projection to SU(N) gauge group

Parameters
----------
maximum_number_of_iteration : int
converegence_criterion : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Projection_Maximum_SU_N*
                  {
                    auto obj = new Projection_Maximum_SU_N;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  py::class_<Projection_Stout_SU3, Projection>(m, "Projection_Stout_SU3", R"pybdoc(
Stout(exponential)-type projection to SU(N) gauge group

Parameters
----------
(none)
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Projection_Stout_SU3*
                  {
                    auto obj = new Projection_Stout_SU3;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    ;

  //! Director_Smear class
  //
  py::class_<Director>(m, "Director", R"pybdoc(
Manager of smeared configurations

Parameters
----------
number_of_smearing : int
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Director*
                  {
                    warn("invalid instantiation of director base class");
                    return nullptr;
                  }), REDIRECT)
    ;

  py::class_<Director_Smear, Director>(m, "Director_Smear")
    .def(py::init([](py::object smear, py::kwargs kwargs) -> Director_Smear*
                  {
                    auto obj = new PyWrap1<Director_Smear, Smear>(smear);

                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("get_Nsmear", &Director_Smear::get_Nsmear, REDIRECT)

    .def("set_config", (void (Director_Smear::*)(Field *)) &Director_Smear::set_config, REDIRECT)
    .def("get_config", [](Director_Smear& self)
         {
           return *(self.get_config());
         }, REDIRECT)
    .def("get_smeared_config", [](Director_Smear& self, int ismear)
         {
           return *(self.getptr_smearedConfig(ismear));
         }, REDIRECT)

    .def("notify_linkv", &Director_Smear::notify_linkv, REDIRECT)
    ;
}
