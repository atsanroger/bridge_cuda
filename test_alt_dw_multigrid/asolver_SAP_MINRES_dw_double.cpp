/*!
        @file    afopr_SAP_MINRES_double.cpp
        @brief   MINRES solver inside a SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2578 $
 */

//====================================================================
#include "asolver_SAP_MINRES_dw.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"

typedef AField<double, QXS>  AField_d;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;

#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"

typedef AField<double, ACCEL>  AField_d;
#endif

#include "lib_alt/Fopr/afopr_dd.h"

#define DEBUG_MINRES

#include "asolver_SAP_MINRES_dw-tmpl.h"

//====================================================================
// explicit instanciation for AField<double>.

#ifdef USE_ALT_QXS
template<>
const std::string ASolver_SAP_MINRES_dw<AField_d>::class_name
                 = "ASolver_SAP_MINRES_dw<AField<double,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string ASolver_SAP_MINRES_dw<AField_d>::class_name
                 = "ASolver_SAP_MINRES_dw<AField<double,ACCEL> >";
#endif

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_SAP_MINRES_dw<AField_d>::register_factory();
}
#endif

template class ASolver_SAP_MINRES_dw<AField_d>;

//============================================================END=====
