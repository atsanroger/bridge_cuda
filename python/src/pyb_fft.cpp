/*!
        @file    pyb_fft.cpp

        @brief   python binding of FFT classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_field.h"
#include "pyb_fft.h"
#include "bridge.h"
#include "Field/field.h"
#include "Tools/fft.h"
#include "Parameters/parameters.h"

void pybind_fft(py::module m)
{
#ifdef USE_FFTWLIB

  py::enum_<FFT::Direction>(m, "FFT_Direction")
    .value("FFT_Forward", FFT::FORWARD)
    .value("FFT_Backward", FFT::BACKWARD)
    .export_values();

  //! FFT class
  //
  py::class_<FFT>(m, "FFT", R"pybdoc(
FFT class

See Also
--------
FFT_xyz_1dim
FFT_xyz_3dim
FFT_3d_local
FFT_3d_parallel1d
FFT_3d_parallel3d
)pybdoc")
    .def(py::init([](const std::string& type_, py::kwargs kwargs) -> FFT*
                  {
                    std::string type = type_;

                    // adjust to factory entries.
                    if (type == "FFT_3d_parallel1d") {
                      type = "FFT_3d_parallel_1dim";
                    } else if (type == "FFT_3d_parallel3d") {
                      type = "FFT_3d_parallel_3dim";
                    }

                    auto obj = FFT::New(type);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }),
         py::arg("type"),
         REDIRECT)

    .def("fft", [](FFT& self, Field& field_out, const Field& field_in, FFT::Direction direction)
         {
           return self.fft(field_out, field_in, direction);
         },
         REDIRECT)
    ;
  
  //! FFT subclasses
  //
  
#define DEFINE_FFT_SUBCLASS(type_)                                      \
  py::class_<type_, FFT>(m, #type_)                                     \
    .def(py::init([](py::kwargs kwargs) -> type_*                       \
                  {                                                     \
                    auto obj = new type_;                               \
                    if (obj) obj->set_parameters(make_parameters(kwargs)); \
                    return obj;                                         \
                  }),                                                   \
         REDIRECT)                                                      \
    ;                                                                   \

  DEFINE_FFT_SUBCLASS(FFT_xyz_1dim)
  DEFINE_FFT_SUBCLASS(FFT_xyz_3dim)

  DEFINE_FFT_SUBCLASS(FFT_3d_local)
#ifdef USE_MPI
  DEFINE_FFT_SUBCLASS(FFT_3d_parallel1d)
  DEFINE_FFT_SUBCLASS(FFT_3d_parallel3d)
#endif

#undef DEFINE_FFT_3d  

#endif
}
