/*!
        @file    asolver_SAP_MINRES_double.cpp
        @brief   MINRES solver inside a SAP solver (Accel version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2225 $
 */


#include "lib_alt/Solver/asolver_SAP_MINRES.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"

//#define DEBUG_MINRES

#include "lib_alt/Solver/asolver_SAP_MINRES-tmpl.h"

//====================================================================
// explicit instanciation for AField<double>.
template<>
const std::string ASolver_SAP_MINRES<AField<double,ACCEL> >::class_name
                              = "ASolver_SAP_MINRES<AField<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_SAP_MINRES<AField<double,ACCEL> >::register_factory();
}
#endif

template class ASolver_SAP_MINRES<AField<double,ACCEL> >;

//============================================================END=====
