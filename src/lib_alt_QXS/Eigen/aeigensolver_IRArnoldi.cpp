/*!
        @file    aeigensolver_IRArnoldi.cpp
        @brief
        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib/Eigen/aeigensolver_IRArnoldi.h"
#include "lib/Eigen/aeigensolver_IRArnoldi-tmpl.h"

#include "lib/Fopr/afopr.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string
AEigensolver_IRArnoldi<AField<double, QXS>,
                       AFopr<AField<double, QXS> > >::class_name
  = "AEigensolver_IRArnoldi<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AField<double, QXS>,
                            AFopr<AField<double, QXS> > >::
               Factory_params::Register("IRArnoldi", create_object);
}
#endif

template class AEigensolver_IRArnoldi<AField<double, QXS>,
                                      AFopr<AField<double, QXS> > >;



// explicit instanciation for AField<float,QXS>.
template<>
const std::string
AEigensolver_IRArnoldi<AField<float, QXS>,
                       AFopr<AField<float, QXS> > >::class_name
  = "AEigensolver_IRArnoldi<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AField<float, QXS>,
                            AFopr<AField<float, QXS> > >::
               Factory_params::Register("IRArnoldi", create_object);
}
#endif

template class AEigensolver_IRArnoldi<AField<float, QXS>,
                                      AFopr<AField<float, QXS> > >;


//============================================================END=====
