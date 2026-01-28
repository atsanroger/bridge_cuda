/*!
        @file    pyb_tools_su_n.cpp

        @brief   python binding of SU(N) matrix and vector classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include <pybind11/operators.h>
#include "pyb_defs.h"
#include "pyb_tools_su_n.h"
#include "bridge.h"
#include "Tools/vec_SU_N.h"
#include "Tools/mat_SU_N.h"

void pybind_tools_su_n(py::module m)
{
  using namespace SU_N;

  py::class_<Vec_SU_N>(m, "Vec_SU_N")
    .def(py::init([](){ return Vec_SU_N(); }), REDIRECT)
    .def(py::init([](int Nc, double r){ return Vec_SU_N(Nc, r); }),
         py::arg("Nc"),
         py::arg("r") = 0.0,
         REDIRECT)

    .def_property_readonly("nc", &Vec_SU_N::nc, REDIRECT)

    .def("set", [](Vec_SU_N& self, const int c, const dcomplex& v) {
        return self.set(c, real(v), imag(v));
      }, REDIRECT)
    .def("get", [](Vec_SU_N& self, const int c) {
        return cmplx(self.r(c), self.i(c));
      }, REDIRECT)

    .def("norm", &Vec_SU_N::norm, REDIRECT)
    .def("dag", &Vec_SU_N::dag, REDIRECT)
    .def("zero", &Vec_SU_N::zero, REDIRECT)
    .def("xI", &Vec_SU_N::xI, REDIRECT)

    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self *= double())
    .def(py::self /= double())
    .def(py::self *= dcomplex())
    .def(py::self /= dcomplex())
    // .def(-py::self)
    ;

  py::class_<Mat_SU_N>(m, "Mat_SU_N")
    .def(py::init<int, double>(),
         py::arg("Nc"),
         py::arg("r") = 0.0,
         REDIRECT)

    .def_property_readonly("nc", &Mat_SU_N::nc, REDIRECT)

    .def("set", [](Mat_SU_N& self, int a, int b, const dcomplex& v) {
        return self.set(a, b, real(v), imag(v));
      },
      py::arg("a"), py::arg("b"), py::arg("v"),
      REDIRECT)
    .def("get", [](Mat_SU_N& self, int a, int b) {
        return cmplx(self.r(a, b), self.i(a, b));
      },
      py::arg("a"), py::arg("b"),
      REDIRECT)

    .def("dag", &Mat_SU_N::dag, REDIRECT)
    .def("ht", &Mat_SU_N::ht, REDIRECT)
    .def("ah", &Mat_SU_N::ah, REDIRECT)
    .def("at", &Mat_SU_N::at, REDIRECT)
    .def("unit", &Mat_SU_N::unit, REDIRECT)
    .def("zero", &Mat_SU_N::zero, REDIRECT)
    .def("xI", &Mat_SU_N::xI, REDIRECT)
    .def("reunit", &Mat_SU_N::reunit, REDIRECT)
    // .def("set_random", &Mat_SU_N::set_random, REDIRECT)

    .def("norm2", &Mat_SU_N::norm2, REDIRECT)

    .def("mult_nn", &Mat_SU_N::mult_nn, REDIRECT)
    .def("mult_nd", &Mat_SU_N::mult_nd, REDIRECT)
    .def("mult_dn", &Mat_SU_N::mult_dn, REDIRECT)
    
    .def("multadd_nn", &Mat_SU_N::multadd_nn, REDIRECT)
    .def("multadd_nd", &Mat_SU_N::multadd_nd, REDIRECT)
    .def("multadd_dn", &Mat_SU_N::multadd_dn, REDIRECT)

    .def("zaxpy", [](Mat_SU_N& self, const dcomplex& a, const Mat_SU_N& x) {
        return self.zaxpy(real(a), imag(a), x);
      }, REDIRECT)
    
    .def(py::self += py::self)
    .def(py::self += double())
    .def(py::self -= py::self)
    .def(py::self -= double())
    .def(py::self *= py::self)
    .def(py::self *= double())
    .def(py::self /= double())
    // .def(-py::self)
    ;

  py::class_<GeneratorSet_Mat_SU_N>(m, "GeneratorSet_Mat_SU_N")
    .def(py::init<int>())
    .def("get_generator", &GeneratorSet_Mat_SU_N::get_generator, REDIRECT)
    ;
}
