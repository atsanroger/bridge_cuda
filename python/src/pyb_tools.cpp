/*!
        @file    pyb_tools.cpp

        @brief   python binding of tools and mathematical functions

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_tools.h"
#include "bridge.h"
#include "Tools/decompose_Hessenberg_Cmplx.h"
#include "Tools/decompose_LUP_Cmplx.h"
#include "Tools/decompose_LU_Cmplx.h"
#include "Tools/decompose_QR_Cmplx.h"
#include "Tools/eigen_QR_Cmplx.h"
#include "Parameters/parameters.h"

#ifdef TRUNK
class PyMath_Sign_Zolotarev : public Math_Sign_Zolotarev
{
 public:
  PyMath_Sign_Zolotarev(int Np, double bmax)
    : Math_Sign_Zolotarev(Np, bmax), Np(Np) {}
 public:
  int Np;
};
#endif

void pybind_module_tools(py::module m)
{
  py::module mm = m.def_submodule("Tools");
  pybind_tools(mm);
}

void pybind_tools(py::module m)
{
  py::class_<Decompose_Hessenberg_Cmplx>(m, "Decompose_Hessenberg_Cmplx")
    .def(py::init<int>(), py::arg("Nin"), REDIRECT)

    .def("set_matrix", &Decompose_Hessenberg_Cmplx::set_matrix, REDIRECT)
    .def("get_Hessenberg", &Decompose_Hessenberg_Cmplx::get_Hessenberg, REDIRECT)
    .def("get_Q", &Decompose_Hessenberg_Cmplx::get_Q, REDIRECT)
    .def("mult_Q", &Decompose_Hessenberg_Cmplx::mult_Q, REDIRECT)
    ;

  py::class_<Decompose_LUP_Cmplx>(m, "Decompose_LUP_Cmplx")
    .def(py::init<int>(), py::arg("N"), REDIRECT)

    .def("set_matrix", &Decompose_LUP_Cmplx::set_matrix, REDIRECT)
    .def("solve", &Decompose_LUP_Cmplx::solve, REDIRECT)
    .def("get_inverse", &Decompose_LUP_Cmplx::get_inverse, REDIRECT)
    .def("mult_inverse", &Decompose_LUP_Cmplx::mult_inverse, REDIRECT)
    .def("determinant", &Decompose_LUP_Cmplx::determinant, REDIRECT)
    ;  

  py::class_<Decompose_LU_Cmplx>(m, "Decompose_LU_Cmplx")
    .def(py::init<int>(), py::arg("N"), REDIRECT)

    .def("set_matrix", &Decompose_LU_Cmplx::set_matrix, REDIRECT)
    .def("solve", &Decompose_LU_Cmplx::solve, REDIRECT)
    .def("get_inverse", &Decompose_LU_Cmplx::get_inverse, REDIRECT)
    .def("mult_inverse", &Decompose_LU_Cmplx::mult_inverse, REDIRECT)
    .def("determinant", &Decompose_LU_Cmplx::determinant, REDIRECT)
    ;  

  py::class_<Decompose_QR_Cmplx>(m, "Decompose_QR_Cmplx")
    .def(py::init<int>(), py::arg("Nin"), REDIRECT)

    .def("set_matrix", &Decompose_QR_Cmplx::set_matrix, REDIRECT)
    .def("get_R", &Decompose_QR_Cmplx::get_R, REDIRECT)
    .def("solve", &Decompose_QR_Cmplx::solve, REDIRECT)
    .def("get_inverse", &Decompose_QR_Cmplx::get_inverse, REDIRECT)
    .def("mult_inverse", &Decompose_QR_Cmplx::mult_inverse, REDIRECT)
    .def("get_Q", &Decompose_QR_Cmplx::get_Q, REDIRECT)
    .def("get_Qu", &Decompose_QR_Cmplx::get_Qu, REDIRECT)
    .def("mult_Q", &Decompose_QR_Cmplx::mult_Q, REDIRECT)
    ;

  py::class_<Eigen_QR_Cmplx>(m, "Eigen_QR_Cmplx")
    .def(py::init<int>(), py::arg("N"), REDIRECT)

    .def("solve", &Eigen_QR_Cmplx::solve, REDIRECT)
    .def("get_R", &Eigen_QR_Cmplx::get_R, REDIRECT)
    .def("get_Q", &Eigen_QR_Cmplx::get_Q, REDIRECT)
    ;


  py::class_<Math_Rational>(m, "Math_Rational")
    .def(py::init([](py::kwargs kwargs)
                  {
                    auto obj = new Math_Rational;
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         REDIRECT)

    .def("set_parameters", [](Math_Rational& self, py::kwargs kwargs)
                  {
                    self.set_parameters(make_parameters(kwargs));
                  },
         REDIRECT)
    .def("get_parameters", [](Math_Rational& self)
         {
           double norm;
           std::vector<double> res;
           std::vector<double> pole;
           self.get_parameters(norm, res, pole);
           return std::make_tuple(norm, res, pole);
         },
         REDIRECT)
    
    .def("func", &Math_Rational::func, REDIRECT)
    ;

#ifdef TRUNK
  py::class_<PyMath_Sign_Zolotarev>(m, "Math_Sign_Zolotarev")
    .def(py::init<int,double>(),
         py::arg("Np"),
         py::arg("bmax"),
         REDIRECT)

    .def("get_sign_parameters", [](PyMath_Sign_Zolotarev& self)
         {
           int np = self.Np;
           std::vector<double> cl(np * 2);
           std::vector<double> bl(np);

           self.get_sign_parameters(cl, bl);

           return std::make_tuple(cl, bl);
         }, REDIRECT)

    .def("sign", &Math_Sign_Zolotarev::sign, REDIRECT)
    ;
#endif

  py::class_<EpsilonTensor>(m, "EpsilonTensor")
    .def(py::init<>())
    .def("index", &EpsilonTensor::epsilon_3_index, REDIRECT)
    .def("value", &EpsilonTensor::epsilon_3_value, REDIRECT)
    ;
}
