/*!
        @file    aeigensolver_IRArnoldi_block.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-04-30 12:32:42 #$
        @version $LastChangedRevision: 2601 $
*/

#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt_Accel/inline/afield_th-inc.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"

#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block-tmpl.h"

// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string
  AEigensolver_IRArnoldi_block<AField<double,ACCEL>,
                         AFopr<AField<double,ACCEL> > >::class_name
  = "AEigensolver_IRArnoldi_block<AField<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AField<double,ACCEL>,
                            AFopr<AField<double, ACCEL> > >::
       Factory_params::Register("IRArnoldi_block", create_object);
}
#endif

template class AEigensolver_IRArnoldi_block<AField<double,ACCEL>,
                                      AFopr<AField<double,ACCEL> > >;

//============================================================END=====
