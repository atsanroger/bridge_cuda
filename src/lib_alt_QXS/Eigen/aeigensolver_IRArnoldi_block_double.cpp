/*!
        @file    aeigensolver_IRArnoldi_block_double.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-04-01 12:13:09 #$
        @version $LastChangedRevision: 2595 $
*/

#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

#include "lib_alt/Fopr/afopr_dd.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"


#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block-tmpl.h"

typedef AField<double, QXS>  AField_d;

template<>
const std::string
AEigensolver_IRArnoldi_block<AField_d, AFopr<AField_d > >::
  class_name = "AEigensolver_IRArnoldi_block<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AField_d>, AFopr<AField_d> >::
  Factory_params::Register("IRArnoldi_block", create_object);
}
#endif

template class AEigensolver_IRArnoldi_block<AField_d, AFopr<AField_d> >;

//============================================================END=====
