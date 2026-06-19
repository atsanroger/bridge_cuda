/*!
        @file    pyb_source.cpp

        @brief   python binding of source classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_source.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Fermion/source.h"
#include "Parameters/parameters.h"

void pybind_source(py::module m)
{
  //! Source
  //
  py::class_<Source>(m, "Source", R"pybdoc(
Sources for linear equation solver
)pybdoc")
    .def(py::init([](const std::string& subtype, py::kwargs kwargs) -> Source*
                  {
                    auto obj = Source::New(subtype);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

    .def("set", [](Source& self,
                   Field& src, int index)
         {
           self.set(src, index);
         },
         py::arg("target"), py::arg("index"),
         REDIRECT
    )
    .def("set", [](Source& self,
                   Field& src, int icolor, int ispin)
         {
           self.set(src, icolor, ispin);
         },
         py::arg("target"), py::arg("icolor"), py::arg("ispin"),
         REDIRECT
    )
    .def("set_all_color", [](Source& self,
                             Field& src, int ispin)
         {
           self.set_all_color(src, ispin);
         },
         py::arg("target"), py::arg("ispin"),
         REDIRECT
    )
    .def("set_all_color_spin", [](Source& self,
                                  Field& src)
         {
           self.set_all_color_spin(src);
         },
         py::arg("target"),
         REDIRECT
    )

    .def("set_parameteres", [](Source& self, py::kwargs kwargs)
         {
           self.set_parameters(make_parameters(kwargs));
         }, REDIRECT
    )
    ;

#define CONSTRUCTOR(T)                                          \
  def(py::init([](py::kwargs kwargs) {                          \
        auto obj = new T;                                       \
        if (obj) obj->set_parameters(make_parameters(kwargs));  \
        return obj;                                             \
      }), REDIRECT)                                             \

  //! Source_Local
  //
  //  parameters:
  //    source_position (int vector)
  //
  py::class_<Source_Local, Source>(m, "Source_Local", R"pybdoc(
Local source for 4-spinor fermion

Parameters
----------
source_position : List [int]
)pybdoc")
    .CONSTRUCTOR(Source_Local)
    .def("set_parameters", [](Source_Local& self, const std::vector<int>& pos)
         {
           self.set_parameters(pos);
         },
         py::arg("source_position"),
         REDIRECT
    )
    ;

  //! Source_Exponential
  //
  //  parameters:
  //    source_position (int vector)
  //    slope (double)
  //    power (double)
  //
  py::class_<Source_Exponential, Source>(m, "Source_Exponential", R"pybdoc(
Exponentially smeared source for 4-spinor fermion

Parameters
----------
source_position : List [int]
slope : float
power : float
)pybdoc")
    .CONSTRUCTOR(Source_Exponential)
    .def("set_parameters", [](Source_Exponential& self,
                              const std::vector<int>& pos,
                              const double slope,
                              const double power)
         {
           self.set_parameters(pos, slope, power);
         },
         py::arg("source_position"),
         py::arg("slope"),
         py::arg("power"),
         REDIRECT
    )
    ;
         
  //! Source_Wall
  //
  //  parameters:
  //    source_position (int vector)
  //
  py::class_<Source_Wall, Source>(m, "Source_Wall", R"pybdoc(
Wall source for 4-spinor fermion

Parameters
----------
source_position : List [int]
)pybdoc")
    .CONSTRUCTOR(Source_Wall)
    .def("set_parameters", [](Source_Wall& self, const std::vector<int>& pos)
         {
           self.set_parameters(pos);
         },
         py::arg("source_position"),
         REDIRECT
    )
    ;
         
  //! Source_MomentumWall
  //
  //  parameters:
  //    source_position (int vector)
  //    source_momentum (int vector)
  //
  py::class_<Source_MomentumWall, Source>(m, "Source_MomentumWall", R"pybdoc(
Momentum wall source for 4-spinor fermion

Parameters
----------
source_position : List [int]
source_momentum : List [int]
)pybdoc")
    .CONSTRUCTOR(Source_MomentumWall)
    .def("set_parameters", [](Source_MomentumWall& self,
                              const std::vector<int>& pos,
                              const std::vector<int>& mom)
         {
           self.set_parameters(pos, mom);
         },
         py::arg("source_position"),
         py::arg("source_momentum"),
         REDIRECT
    )
    ;
         
  //! Source_Random
  //
  //  parameters:
  //    source_position (int vector)
  //    source_momentum (int vector)
  //    noise_type (string) = { Gaussian, U1, Z2 }
  //
  //  extra method:
  //    set_all_space_time(src, icolor, ispin = -1)
  //      fill all spin index when ispin = -1
  //
  py::class_<Source_Random, Source>(m, "Source_Random", R"pybdoc(
Random noise source at a given time slice

Parameters
----------
source_position : List [int]
source_momentum : List [int]
noise_type : string
)pybdoc")
    .CONSTRUCTOR(Source_Random)
    .def("set_parameters", [](Source_Random& self,
                              const std::vector<int>& pos,
                              const std::vector<int>& mom,
                              const std::string& noise)
         {
           self.set_parameters(pos, mom, noise);
         },
         py::arg("source_position"),
         py::arg("source_momentum"),
         py::arg("noise_type"),
         REDIRECT
    )

    .def("set_all_space_time", [](Source_Random& self,
                                  Field& src, int icolor, int ispin)
         {
           if (ispin == -1) {
             self.set_all_space_time(src, icolor);
           } else {
             self.set_all_space_time(src, icolor, ispin);
           }
         },
         py::arg("target"), py::arg("icolor"), py::arg("ispin") = -1,
         REDIRECT
    )
    ;

#ifdef TRUNK
  //! Source_Staggered_Wall
  //
  py::class_<Source_Staggered_Wall>(m, "Source_Staggered_Wall", R"pybdoc(
Wall source for staggered fermion

Parameters
----------
source_position : List [int]
)pybdoc")
    .CONSTRUCTOR(Source_Staggered_Wall)
    .def("set_parameters", [](Source_Staggered_Wall& self,
                              const int pos)
         {
           self.set_parameters(pos);
         },
         py::arg("source_position"),
         REDIRECT)

    .def("set", &Source_Staggered_Wall::set,
         py::arg("src"),
         py::arg("icolor"),
         py::arg("isource"),
         REDIRECT)
    ;
#endif

#undef CONSTRUCTOR
}
