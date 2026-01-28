/*!
        @file    pyb_eigen.cpp

        @brief   python binding of Eigensolver classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <tuple>
#include <vector>
#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_eigen.h"
#include "bridge.h"
#include "Field/field.h"
#include "Eigen/eigensolver.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

void pybind_eigen(py::module m)
{
  //! Eivensolver
  //
  py::class_<Eigensolver>(m, "Eigensolver", R"pybdoc(
Eigenvalue solver class

See Also
--------
Eigensolver_IRLanczos
)pybdoc")
    .def(py::init([](){
          warn("invalid instantiation of eigensolver base class");
          return nullptr;
        }), REDIRECT)
    ;

  //! IR Lanczos algorithm
  //
  py::class_<Eigensolver_IRLanczos, Eigensolver>(m, "Eigensolver_IRLanczos",
      R"pybdoc(
Eigenvalue solver with Implicitly Restarted Lanczos algorithm

Parameters
----------
eigensolver_mode : string
number_of_wanted_eigenvectors : int
number_of_working_eigenvectors : int
maximum_number_of_iteration : int
convergence_criterion_squared : float
threshold_value : float
)pybdoc")
    .def(py::init([](py::object fopr, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<Eigensolver_IRLanczos, Fopr>(fopr);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("solve", [](Eigensolver_IRLanczos& self, const Field& b, const int work_size)
         {
           int nin = b.nin();
           int nvol = b.nvol();
           int nex = b.nex();

           std::vector<Field>  vk(work_size, Field(nin,nvol,nex));
           std::vector<double> tda(work_size);

           int nsbt = -1;
           int nconv = -1;

           self.solve(tda, vk, nsbt, nconv, b);

           return std::make_tuple(vk, tda, nsbt, nconv);
         },
         py::arg("source"),
         py::arg("work_size")=100,
         R"pybdoc(
Arguments
---------
source : Field
    source vector
work_size : int, optional
    number of work vectors

Return values
-------------
tuple of the following items
    list of eigenvectors : list [Field]
    list of eigenvalues : list [float]
    number of eigenpairs : int
    number of iterations : int
)pybdoc",
         REDIRECT
    )

    .def("set_parameters", [](Eigensolver_IRLanczos& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT)
    ;
}
