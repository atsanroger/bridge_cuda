/*!
        @file    pyb_qsuscept.cpp

        @brief   python binding of quark number susceptibility measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_qsuscept.h"
#include "bridge.h"
#include "Measurements/Fermion/quarkNumberSusceptibility_Wilson.h"
#include "Measurements/Fermion/noiseVector.h"
#include "Measurements/Fermion/noiseVector_Z2.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

void pybind_measure_qsuscept(py::module m)
{
  //! QuarkNumberSusceptibility
  //
  py::class_<QuarkNumberSusceptibility_Wilson>(m, "QuarkNumberSusceptibility_Wilson", R"pybdoc(
Quark number susceptibility for the Wilson-type fermion

Arguments
---------
fopr
fprop
noise_vector

Parameters
----------
number_of_noises : int
)pybdoc")
    .def(py::init([](py::object fopr,
                     py::object fprop,
                     py::object noise,
                     py::kwargs kwargs)
                  {
                    auto obj = new PyWrap3<QuarkNumberSusceptibility_Wilson,
                                           Fopr,
                                           Fprop,
                                           NoiseVector>(fopr, fprop, noise);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("op"),
         py::arg("fprop"),
         py::arg("noise_vector"),
         REDIRECT)
    .def("set_parameters", [](QuarkNumberSusceptibility_Wilson& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         },
         REDIRECT)
    .def("measure", &QuarkNumberSusceptibility_Wilson::measure)
    ;

  //! NoiseVector
  //
  py::class_<NoiseVector>(m, "NoiseVector", R"pybdoc(
Noise vector generator
)pybdoc")
    .def(py::init([]() -> NoiseVector*
                  {
                    warn("invalid instantiation of NoiseVector base class");
                    return nullptr;
                  }), REDIRECT)

    .def("set", &NoiseVector::set, REDIRECT)
    ;

  //! Z2 NoiseVector
  //
  py::class_<NoiseVector_Z2, NoiseVector>(m, "NoiseVector_Z2")
    .def(py::init([](py::object rand, py::kwargs kwargs)
                  {
                    auto obj = new PyWrap1<NoiseVector_Z2, RandomNumbers>(rand);
                    // if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("rand"),
         REDIRECT)
    ;
          
}
