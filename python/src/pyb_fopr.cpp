/*!
        @file    pyb_fopr.cpp

        @brief   python binding of Fermion operator classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_fopr.h"
#include "pyb_smear.h"
#include "bridge.h"
#include "Field/field.h"
#include "Fopr/fopr.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

// create concrete object of Wilson type
//
Fopr *make_fopr(const std::string& subtype, py::dict args)
{
  if (subtype == "Wilson"
      || subtype == "Wilson_Isochemical"
      || subtype == "WilsonGeneral"
      || subtype == "Clover"
      || subtype == "Clover_Isochemical"
      || subtype == "CloverGeneral"
      || subtype == "Wilson_eo"
      || subtype == "Clover_eo"
#ifdef TRUNK
      || subtype == "Wilson_TwistedMass"
#endif
  ) {

    std::string repr = args.contains("gamma_matrix_type")
      ? args["gamma_matrix_type"].cast<std::string>()
      : "Dirac";

    auto obj = Fopr::New(subtype, repr);
    if (obj) obj->set_parameters(make_parameters(args));

    return obj;
  }

  warn("Fopr: unsupported fermion type.");
  return nullptr;
}

void pybind_fopr(py::module m)
{

  //! Fopr base class
  //    defines interface to Fopr subtypes
  //
  py::class_<Fopr>(m, "Fopr", R"pybdoc(
Fermion Operator class

Returns an instance of generic Fermion operator class. The particular
type is specified by `type` argument.

See Also
--------
Fopr_Wilson
Fopr_Wilson_eo
Fopr_Wilson_Isochemical
Fopr_WilsonGeneral
Fopr_Clover
Fopr_Clover_eo
Fopr_Clover_Isochemical
Fopr_CloverGeneral
)pybdoc"
#ifdef TRUNK
R"pybdoc(
Fopr_Wilson_TwistedMass
Fopr_Staggered_eo
Fopr_Domainwall
Fopr_Overlap
Fopr_Overlap_5d
)pybdoc"
#endif
R"pybdoc(
Fopr_Rational
Fopr_Chebyshev
Fopr_Smeared
Fopr_Smeared_eo
Fopr_CRS
)pybdoc"
  )
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) {
          return make_fopr(subtype, kwargs);
        }),
      py::arg("type"), REDIRECT)

    .def("mult", [](Fopr& self, Field& w, const Field& v) -> void
         {
           return self.mult(w, v);
         }, REDIRECT)
    .def("mult_dag", [](Fopr& self, Field& w, const Field& v) -> void
         {
           return self.mult_dag(w, v);
         }, REDIRECT)

    .def("set_mode", &Fopr::set_mode, REDIRECT)
    .def("get_mode", &Fopr::get_mode, REDIRECT)

    .def("set_config", (void(Fopr::*)(Field*)) &Fopr::set_config, REDIRECT)
    
    .def("field_nin", &Fopr::field_nin, REDIRECT)
    .def("field_nvol", &Fopr::field_nvol, REDIRECT)
    .def("field_nex", &Fopr::field_nex, REDIRECT)

    .def("flop_count", (double (Fopr::*)()) &Fopr::flop_count, REDIRECT)
    ;

#define DEFINE_WILSON_TYPE_CLASS(type_, name_, docstring_)              \
  py::class_<type_, Fopr>(m, #type_, docstring_)                        \
    .def(py::init([](py::kwargs kwargs) {                               \
          std::string repr = kwargs.contains("gamma_matrix_type")       \
            ? kwargs["gamma_matrix_type"].cast<std::string>()           \
            : "Dirac";                                                  \
          auto obj = new type_(repr);                                   \
          if (obj) obj->set_parameters(make_parameters(kwargs));        \
          return obj;                                                   \
        }), REDIRECT)                                                   \
    .def("mult_gm5", [](type_& self, Field& w, const Field& v) -> void  \
         {                                                              \
           return self.mult_gm5(w, v);                                  \
         }, REDIRECT)                                                   \
    ;

  DEFINE_WILSON_TYPE_CLASS(Fopr_Wilson, Wilson, R"pybdoc(
Wilson fermion operator

Parameters
----------
hopping_parameter : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")

  DEFINE_WILSON_TYPE_CLASS(Fopr_Clover, Clover, R"pybdoc(
Clover fermion operator

Parameters
----------
hopping_parameter : float
clover_coefficient : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")

  DEFINE_WILSON_TYPE_CLASS(Fopr_Wilson_Isochemical, Wilson_Isochemical, R"pybdoc(
Wilson fermion operator with isospin chemical potential

Parameters
----------
hopping_parameter : float
isospin_chemical_potential : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")

  DEFINE_WILSON_TYPE_CLASS(Fopr_Clover_Isochemical, Clover_Isochemical, R"pybdoc(
Clover fermion operator with isospin chemical potential

Parameters
----------
hopping_parameter : float
clover_coefficient : float
isospin_chemical_potential : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")

  DEFINE_WILSON_TYPE_CLASS(Fopr_WilsonGeneral, WilsonGeneral, R"pybdoc(
Wilson general fermion operator

Parameters
----------
hopping_parameter_spatial : float
hopping_parameter_temporal : float
dispersion_parameter_spatial : float
Wilson_parameter_spatial : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")

  DEFINE_WILSON_TYPE_CLASS(Fopr_CloverGeneral, CloverGeneral, R"pybdoc(
Clover general fermion operator

Parameters
----------
hopping_parameter_spatial : float
hopping_parameter_temporal : float
dispersion_parameter_spatial : float
Wilson_parameter_spatial : float
clover_coefficient_spatial : float
clover_coefficient_temporal : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")

#ifdef TRUNK
  DEFINE_WILSON_TYPE_CLASS(Fopr_Wilson_TwistedMass, Wilson_TwistedMass, R"pybdoc(
Twisted-mass Wilson fermion operator

Parameters
----------
hopping_parameter : float
twisted_mass : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")
#endif

#undef DEFINE_WILSON_TYPE_CLASS

  //! Composite classes

  //! Fopr_Rational for rational approximation
  //
  py::class_<Fopr_Rational, Fopr>(m, "Fopr_Rational", R"pybdoc(
Fermion operator for rational approximation

Parameters
----------
number_of_poles : int
exponent_numerator : int
exponent_denomicator : int
lower_bound : float
upper_bound : float
maximum_number_of_iteration : int
convergence_criterion_squared : float
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) {
          auto obj = new PyWrap1<Fopr_Rational, Fopr>(arg);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;

  //! Fopr_Chebyshev for Chebyshev approximation
  //
  py::class_<Fopr_Chebyshev, Fopr>(m, "Fopr_Chebyshev", R"pybdoc(
Chebyshev polynomial of fermion operator

Parameters
----------
degree_of_polynomial : int
threshold_value : float
upper_bound : float
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) {
          auto obj = new PyWrap1<Fopr_Chebyshev, Fopr>(arg);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    .def("mult", [](Fopr_Chebyshev& self, double x) -> double
         {
           return self.mult(x);
         }, REDIRECT)
    ;

  //! Fopr_Smeared for decorating with smearing
  //
  py::class_<Fopr_Smeared, Fopr>(m, "Fopr_Smeared", R"pybdoc(
Smeared fermion operator

Parameters
----------
(none)
)pybdoc")
    .def(py::init([](py::object fopr, py::object smear, py::kwargs kwargs) {
          auto obj = new PyWrap2<Fopr_Smeared, Fopr, Director_Smear>(fopr, smear);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;

  //! Fopr_CRS for handling Dirac operator in CRS data format
  //
  py::class_<Fopr_CRS, Fopr>(m, "Fopr_CRS", R"pybdoc(
Fermion operator with CRS matrix format

Parameters
----------
(none)
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) -> Fopr_CRS* {

          const char *tp = arg.ptr()->ob_type->tp_name;

          if (strcmp(tp, "str") == 0) {
            // string passed as argument: Fopr_CRS(filename)

            auto obj = new Fopr_CRS(arg.cast<std::string>());
            if (obj) obj->set_parameters(make_parameters(kwargs));

            return obj;

          } else {
            // assume Fopr passed as argument: Fopr_CRS(Fopr*)

            auto obj = new PyWrap1<Fopr_CRS, Fopr>(arg);
            if (obj) obj->set_parameters(make_parameters(kwargs));

            return obj;
          }
        }), REDIRECT)

    .def("write_matrix", &Fopr_CRS::write_matrix, REDIRECT)
    ;


  //! Fopr_eo - base class for even-odd preconditioned operators
  //
  py::class_<Fopr_eo, Fopr>(m, "Fopr_eo")
    .def(py::init([](py::kwargs kwargs)
                  {
                    warn("invalid instantiation of Fopr_eo base class");
                    return nullptr;
                  }),
         REDIRECT)
    ;

  //! Fopr_Wilson_eo
  //
  py::class_<Fopr_Wilson_eo, Fopr_eo>(m, "Fopr_Wilson_eo", R"pybdoc(
Even-odd Wilson fermion operator

Parameters
----------
hopping_parameter : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          std::string repr = kwargs.contains("gamma_matrix_type")
            ? kwargs["gamma_matrix_type"].cast<std::string>()
            : "Dirac";
          auto obj = new Fopr_Wilson_eo(repr);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    .def("mult_gm5", [](Fopr_Wilson_eo& self, Field& w, const Field& v) -> void
         {
           return self.mult_gm5(w, v);
         }, REDIRECT)
    ;          

  //! Fopr_Clover_eo
  //
  py::class_<Fopr_Clover_eo, Fopr_eo>(m, "Fopr_Clover_eo", R"pybdoc(
Even-odd Clover fermion operator

Parameters
----------
hopping_parameter : float
clover_coefficient : float
boundary_condition : List [int]
gamma_matrix_type : string, optional (default='Dirac')
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          std::string repr = kwargs.contains("gamma_matrix_type")
            ? kwargs["gamma_matrix_type"].cast<std::string>()
            : "Dirac";
          auto obj = new Fopr_Clover_eo(repr);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    .def("mult_gm5", [](Fopr_Clover_eo& self, Field& w, const Field& v) -> void
         {
           return self.mult_gm5(w, v);
         }, REDIRECT)
    ;          

  //! Fopr_Smeared_eo for decorating with smearing
  //
  py::class_<Fopr_Smeared_eo, Fopr_eo>(m, "Fopr_Smeared_eo", R"pybdoc(
Smeared fermion operator with even-odd preconditioning

Parameters
----------
(none)
)pybdoc")
    .def(py::init([](py::object fopr, py::object smear, py::kwargs kwargs) {
          auto obj = new PyWrap2<Fopr_Smeared_eo, Fopr_eo, Director_Smear>(fopr, smear);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)
    ;

#ifdef TRUNK
  py::class_<Fopr_Staggered_eo, Fopr>(m, "Fopr_Staggered_eo", R"pybdoc(
Even-odd staggered fermion operator

Parameters
----------
quark_mass : float
boundary_condition : List [int]
)pybdoc")
    .def(py::init([](py::kwargs kwargs) {
          auto obj = new Fopr_Staggered_eo;
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("set_config", [](Fopr_Staggered_eo& self, Field_G* ueo)
         {
           return self.set_config(ueo);
         }, REDIRECT)
    ;

  py::class_<Fopr_Domainwall, Fopr>(m, "Fopr_Domainwall", R"pybdoc(
Domain-wall fermion operator

Parameters
----------
quark_mass : float
domain_wall_height : float
extent_of_5th_dimension : int
boundary_condition : List [int]
gamma_matrix_type : string, optional
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) {
          auto obj = new PyWrap1<Fopr_Domainwall, Fopr_Wilson>(arg);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("mult_chiral_projection_4d", &Fopr_Domainwall::mult_chproj_4d, REDIRECT)
    ;

  py::class_<Fopr_Overlap, Fopr>(m, "Fopr_Overlap", R"pybdoc(
Overlap fermion operator

Parameters
----------
quark_mass : float
domain_wall_height : float
number_of_poles : int
lower_bound : float
upper_bound : float
maximum_number_of_iteration : int
convergence_criterion_squared : float
boundary_condition : List [int]
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) {
          auto obj = new PyWrap1<Fopr_Overlap, Fopr_Wilson>(arg);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("set_lowmodes", &Fopr_Overlap::set_lowmodes, REDIRECT)
    ;

  py::class_<Fopr_Overlap_5d, Fopr>(m, "Fopr_Overlap_5d", R"pybdoc(
Overlap fermion operator in 5-dimensional formulation

Parameters
----------
quark_mass : float
domain_wall_height : float
number_of_poles : int
lower_bound : float
upper_bound : float
maximum_number_of_iteration : int
convergence_criterion_squared : float
boundary_condition : List [int]
)pybdoc")
    .def(py::init([](py::object arg, py::kwargs kwargs) {
          auto obj = new PyWrap1<Fopr_Overlap_5d, Fopr_Wilson_eo>(arg);
          if (obj) obj->set_parameters(make_parameters(kwargs));
          return obj;
        }), REDIRECT)

    .def("set_lowmodes", &Fopr_Overlap_5d::set_lowmodes, REDIRECT)
    ;

#endif

}
