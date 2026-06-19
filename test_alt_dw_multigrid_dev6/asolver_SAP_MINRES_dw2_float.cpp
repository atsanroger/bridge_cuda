/*!
        @file    afopr_SAP_MINRES_dw2_float.cpp
        @brief   MINRES solver inside a SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2604 $
 */

//====================================================================
#include "asolver_SAP_MINRES_dw2.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt/Fopr/afopr_dd.h"

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"
typedef AField<float, QXS>  AField_f;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;

#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"

typedef AField<float, ACCEL>  AField_f;
#endif

#define DEBUG_MINRES

#include "asolver_SAP_MINRES_dw2-tmpl.h"

//====================================================================
// explicit instanciation for AField<float>.
#ifdef USE_ALT_QXS
template<>
const std::string ASolver_SAP_MINRES_dw2<AField_f>::class_name
                  = "ASolver_SAP_MINRES_dw2<AField<float,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string ASolver_SAP_MINRES_dw2<AField_f>::class_name
                  = "ASolver_SAP_MINRES_dw2<AField<float,ACCEL> >";
#endif

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_SAP_MINRES_dw2<AField_f>::register_factory();
}
#endif

template class ASolver_SAP_MINRES_dw2<AField_f>;

//============================================================END=====
