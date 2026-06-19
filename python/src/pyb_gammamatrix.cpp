/*!
        @file    pyb_gammamatrix.cpp

        @brief   python binding of gamma matrices

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_params.h"
#include "pyb_gammamatrix.h"
#include "bridge.h"
#include "Tools/gammaMatrix.h"
#include "Tools/gammaMatrixSet.h"
#include "Parameters/parameters.h"

void pybind_gammamatrix(py::module m)
{
  //! GammaMatrixSet::GMspecies
  //
  py::enum_<GammaMatrixSet::GMspecies>(m, "GMspecies")

#define DEFINE_ITEM(x) value(#x, GammaMatrixSet:: x)

    .DEFINE_ITEM(UNITY)

    .DEFINE_ITEM(GAMMA1)
    .DEFINE_ITEM(GAMMA2)
    .DEFINE_ITEM(GAMMA3)
    .DEFINE_ITEM(GAMMA4)

    .DEFINE_ITEM(GAMMA5)

    .DEFINE_ITEM(GAMMA51)
    .DEFINE_ITEM(GAMMA52)
    .DEFINE_ITEM(GAMMA53)
    .DEFINE_ITEM(GAMMA54)

    .DEFINE_ITEM(GAMMA15)
    .DEFINE_ITEM(GAMMA25)
    .DEFINE_ITEM(GAMMA35)
    .DEFINE_ITEM(GAMMA45)

    .DEFINE_ITEM(SIGMA12)
    .DEFINE_ITEM(SIGMA23)
    .DEFINE_ITEM(SIGMA31)

    .DEFINE_ITEM(SIGMA41)
    .DEFINE_ITEM(SIGMA42)
    .DEFINE_ITEM(SIGMA43)

    .DEFINE_ITEM(CHARGECONJG)

#undef DEFINE_ITEM

    .export_values()
    ;


  //! GammaMatrix
  //
  //  parameters:
  //
  py::class_<GammaMatrix>(m, "GammaMatrix")
    .def(py::init<>())

#define DEFINE_METHOD(x) def(#x, &GammaMatrix:: x)

    .DEFINE_METHOD(set)
    .DEFINE_METHOD(set_values)

    .def("mult",   (GammaMatrix (GammaMatrix::*)(GammaMatrix) const) &GammaMatrix::mult)
    .def("mult_i", (GammaMatrix (GammaMatrix::*)(GammaMatrix) const) &GammaMatrix::mult_i)

    .def("mult", (GammaMatrix (GammaMatrix::*)(int) const) &GammaMatrix::mult)

    .DEFINE_METHOD(value)
    .DEFINE_METHOD(index_c)
    .DEFINE_METHOD(value_r)
    .DEFINE_METHOD(value_i)

    .DEFINE_METHOD(print)

#undef DEFINE_METHOD
    ;


  //! GammaMatrixSet
  //
  //  parameters:
  //
  py::class_<GammaMatrixSet>(m, "GammaMatrixSet")
    .def(py::init([](const std::string& subtype) -> GammaMatrixSet*
                  {
                    auto obj = GammaMatrixSet::New(subtype);
                    return obj;
                  }))

#define DEFINE_METHOD(x) def(#x, &GammaMatrixSet:: x)

    // .DEFINE_METHOD(init_GM)
    .DEFINE_METHOD(get_GM)

    .DEFINE_METHOD(print)

#undef DEFINE_METHOD
    ;

}
