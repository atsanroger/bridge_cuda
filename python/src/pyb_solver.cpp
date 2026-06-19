/*!
        @file    pyb_solver.cpp

        @brief   python binding of linear equation solver classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <tuple>
#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_fopr.h"
#include "pyb_solver.h"
#include "bridge.h"
#include "Field/field.h"
#include "Fopr/fopr.h"
#include "Solver/solver.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

Solver* generate_solver(const std::string& subtype, py::object fopr, py::kwargs kwargs)
{
  Solver *obj = nullptr;

  // hand-written factory
  if (subtype == "CG") {
    obj = new PyWrap1<Solver_CG, Fopr>(fopr);
  } else if (subtype == "CGNE") {
    obj = new PyWrap1<Solver_CGNE, Fopr>(fopr);
  } else if (subtype == "CGNR") {
    obj = new PyWrap1<Solver_CGNR, Fopr>(fopr);
  } else if (subtype == "BiCGStab_Cmplx") {
    obj = new PyWrap1<Solver_BiCGStab_Cmplx, Fopr>(fopr);
  } else if (subtype == "BiCGStab_L_Cmplx") {
    obj = new PyWrap1<Solver_BiCGStab_L_Cmplx, Fopr>(fopr);
  } else if (subtype == "BiCGStab_DS_L_Cmplx") {
    obj = new PyWrap1<Solver_BiCGStab_DS_L_Cmplx, Fopr>(fopr);
  } else if (subtype == "BiCGStab_IDS_L_Cmplx") {
    obj = new PyWrap1<Solver_BiCGStab_IDS_L_Cmplx, Fopr>(fopr);
  } else if (subtype == "GMRES_m_Cmplx") {
    obj = new PyWrap1<Solver_GMRES_m_Cmplx, Fopr>(fopr);
  }

  if (obj) obj->set_parameters(make_parameters(kwargs));

  return obj;
}

void pybind_solver(py::module m)
{
  //! Solver classes
  //
  py::class_<Solver>(m, "Solver", R"pybdoc(
Solver(type, op, params)

Returns an instance of generic Solver class. The particular algorithm 
is specified by `type` argument. `op` refers to a fermion operator.

See Also
--------
Solver_CG
Solver_CGNE
Solver_CGNR
Solver_BiCGStab_Cmplx
Solver_BiCGStab_L_Cmplx
Solver_BiCGStab_DS_L_Cmplx
Solver_BiCGStab_IDS_L_Cmplx
Solver_GMRES_m_Cmplx

)pybdoc")
    .def(py::init([](const std::string& subtype, py::object arg, py::kwargs kwargs) -> Solver*
                  {
                    return generate_solver(subtype, arg, kwargs);
                  }
         ), py::arg("type"), py::arg("op"),
         REDIRECT
    )
    .def("solve", [](Solver& self, Field& sol, const Field& src) -> std::tuple<int,double>
         {
           int nconv = 0;
           double diff = 0.0;
           self.solve(sol, src, nconv, diff);
           return std::tuple<int,double>(nconv, diff);
         },
         R"pybdoc(
solve(solution, source) -> (nconv, diff)

Find a solution to the linear equation for the source, and returns 
a tuple containing the number of iterations and the residue.
)pybdoc",
         REDIRECT)
    ;

#define DEFINE_CLASS(T,P)                                               \
  py::class_<T, Solver>(m, #T, P)                                       \
    .def(py::init([](py::object arg, py::kwargs kwargs)                 \
                  {                                                     \
                    auto obj = new PyWrap1<T, Fopr>(arg);               \
                    if (obj) obj->set_parameters(make_parameters(kwargs)); \
                    return obj;                                         \
                  }), py::arg("op"), REDIRECT)

  DEFINE_CLASS(Solver_CG, R"pybdoc(
Conjugate Gradient (CG) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
)pybdoc");

  DEFINE_CLASS(Solver_CGNE, R"pybdoc(
Conjugate Gradient Normal Equation (CGNE) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
)pybdoc");

  DEFINE_CLASS(Solver_CGNR, R"pybdoc(
Conjugate Gradient Normal Residual (CGNR) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
)pybdoc");

  DEFINE_CLASS(Solver_BiCGStab_Cmplx, R"pybdoc(
BiCGStab solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
Omega_tolerance : float
use_initial_guess : {'true', 'false'}
)pybdoc");

  DEFINE_CLASS(Solver_BiCGStab_L_Cmplx, R"pybdoc(
BiCGStab(L) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
Omega_tolerance : float
number_of_orthonormal_vectors : int
)pybdoc");

  DEFINE_CLASS(Solver_BiCGStab_DS_L_Cmplx, R"pybdoc(
BiCGStab(DS_L) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
Omega_tolerance : float
number_of_orthonormal_vectors : int
tolerance_for_DynamicSelection_of_L : float
)pybdoc");

  DEFINE_CLASS(Solver_BiCGStab_IDS_L_Cmplx, R"pybdoc(
BiCGStab(IDS_L) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
Omega_tolerance : float
number_of_orthonormal_vectors : int
tolerance_for_DynamicSelection_of_L : float
)pybdoc");

  DEFINE_CLASS(Solver_GMRES_m_Cmplx, R"pybdoc(
GMRES(m) solver

Parameters
----------
maximum_number_of_iteration : int
maximum_number_of_restart : int
convergence_criterion_squared : float
use_initial_guess : {'true', 'false'}
number_of_orthonormal_vectors : int
)pybdoc");

#undef DEFINE_CLASS

  //! Shiftsolver class
  //
  py::class_<Shiftsolver>(m, "Shiftsolver", R"pybdoc(
Parameters
----------
maximum_number_of_iteration : int
convergence_criterion_squared : float
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Shiftsolver*
                  {
                    warn("invalid instantiation of shift solver base class");
                    return nullptr;
                  }),
         REDIRECT)

    .def("solve", &Shiftsolver::solve, REDIRECT)

    .def("set_parameters", [](Shiftsolver& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         },
         REDIRECT)
    ;

  py::class_<Shiftsolver_CG, Shiftsolver>(m, "Shiftsolver_CG")
    .def(py::init([](py::object fopr, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<Shiftsolver_CG, Fopr>(fopr);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("op"),
         REDIRECT)
    ;

}
