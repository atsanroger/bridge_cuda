/*!
        @file    pyb_correlator.cpp

        @brief   python binding of correlator measurement classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_correlator.h"
#include "bridge.h"
#include "Field/field.h"
#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/corr4pt_4spinor.h"
#include "Tools/gammaMatrixSet.h"
#include "Parameters/parameters.h"
#include "pyb_wrap.h"

void pybind_correlator(py::module m)
{

  //! 2pt Correlator
  //
  //  parameters:
  //    filename_output (string)
  //

  py::class_<Corr2pt_4spinor>(m, "Corr2pt_4spinor", R"pybdoc(
Two-point correlator for Wiilson-type fermions

Parameters
----------
filename_output : string
)pybdoc")
    .def(py::init([](py::object gmset, py::kwargs kwargs) -> Corr2pt_4spinor*
                  {
                    auto obj = new PyWrap1<Corr2pt_4spinor, GammaMatrixSet>(gmset);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)
    
#define DEFINE_METHOD(x) def(#x, &Corr2pt_4spinor:: x, REDIRECT)
    .DEFINE_METHOD(meson_all)
    .DEFINE_METHOD(meson_correlator)
    .DEFINE_METHOD(meson_momentum_all)
    .DEFINE_METHOD(meson_momentum_correlator)

    .DEFINE_METHOD(meson_correlator_x)
    .DEFINE_METHOD(meson_momentum_correlator_x)

    .DEFINE_METHOD(proton_test)
    .DEFINE_METHOD(proton_correlator)

    .DEFINE_METHOD(proton_correlator_x)
#undef DEFINE_METHOD
    ;


  //! 4pt Correlator
  //
  //  parameters:
  //    filename_output (string)
  //

  py::class_<Corr4pt_4spinor>(m, "Corr4pt_4spinor", R"pybdoc(
Four-point correlator for Wiilson-type fermions

Parameters
----------
filename_output : string
)pybdoc")
    .def(py::init([](py::object gmset, py::kwargs kwargs) -> Corr4pt_4spinor*
                  {
                    auto obj = new PyWrap1<Corr4pt_4spinor, GammaMatrixSet>(gmset);
                    if (obj) obj->set_parameters(make_parameters(kwargs));
                    return obj;
                  }), REDIRECT)

#define DEFINE_METHOD(x) def(#x, &Corr4pt_4spinor:: x, REDIRECT)
    .DEFINE_METHOD(meson_all)
    .DEFINE_METHOD(meson_correlator)
    .DEFINE_METHOD(meson_momentum_all)
    .DEFINE_METHOD(meson_momentum_correlator)
#undef DEFINE_METHOD
    ;


#ifdef TRUNK

  //! 2pt Correlator for Staggered Fermion
  //
  py::class_<Corr2pt_Staggered>(m, "Corr2pt_Staggered", R"pybdoc(
Two-point correlator for staggered fermions

Parameters
----------
(none)
)pybdoc")
    .def(py::init([](py::kwargs kwargs) -> Corr2pt_Staggered*
                  {
                    auto obj = new Corr2pt_Staggered;
                    if (obj) {
                      Bridge::VerboseLevel vl = kwargs.contains("verbose_level")
                        ? Bridge::BridgeIO::set_verbose_level(
                          kwargs["verbose_level"].cast<std::string>())
                        : CommonParameters::Vlevel();

                      obj->set_parameter_verboselevel(vl);
                    }
                    return obj;
                  }),
         REDIRECT)

    .def("meson", [](Corr2pt_Staggered& self,
                     std::vector<Field_F_1spinor>& sq1,
                     std::vector<Field_F_1spinor>& sq2)
         {
           std::vector<double> corr;
           self.meson(corr, sq1, sq1);
           return corr;
         }, REDIRECT)
    ;

#endif
}
