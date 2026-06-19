/*!
        @file    aeigensolver_IRLanczos.cpp
        @brief
        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Eigen/aeigensolver_IRLanczos.h"
#include "lib/Eigen/aeigensolver_IRLanczos-tmpl.h"

#include "lib/Fopr/afopr.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

typedef AField<double,SIMD2> AFIELD;
typedef AFopr<AField<double,SIMD2> > AFOPR;

// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string
  AEigensolver_IRLanczos<AField<double,SIMD2>,
                         AFopr<AField<double,SIMD2> > >::class_name
  = "AEigensolver_IRLanczos<AField<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AFIELD,AFOPR>::
    Factory_params::Register("IRLanczos", create_object);
}
#endif


template class AEigensolver_IRLanczos<AField<double,SIMD2>,
                                      AFopr<AField<double,SIMD2> > >;

//============================================================END=====
