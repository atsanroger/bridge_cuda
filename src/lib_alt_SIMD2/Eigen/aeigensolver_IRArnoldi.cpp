/*!
        @file    aeigensolver_IRArnoldi.cpp
        @brief
        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Eigen/aeigensolver_IRArnoldi.h"
#include "lib/Eigen/aeigensolver_IRArnoldi-tmpl.h"

#include "lib/Fopr/afopr.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string
  AEigensolver_IRArnoldi<AField<double,SIMD2>,
                         AFopr<AField<double,SIMD2> > >::class_name
  = "AEigensolver_IRArnoldi<AField<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AField<double,SIMD2>,
                            AFopr<AField<double, SIMD2> > >::
              Factory_params::Register("IRArnoldi", create_object);
}
#endif

template class AEigensolver_IRArnoldi<AField<double,SIMD2>,
                                      AFopr<AField<double,SIMD2> > >;

//============================================================END=====
