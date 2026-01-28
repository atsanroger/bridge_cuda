/*!

        @file    $Id: MultiGrid_Domainwall_float.cpp #$

        @brief   MultiGrid operation for Domainwall fermion

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2595 $

 */

//====================================================================

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"

typedef AField<float, QXS> AField_f;
#endif

#ifdef USE_ALT_ACCEL
typedef float real_t;

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"

typedef AField<float, ACCEL> AField_f;
#endif

// template for MultiGrid_Clover
//#include "lib_alt/Solver/MultiGrid_Domainwall.h"
//#include "lib_alt/Solver/MultiGrid_Domainwall-tmpl.h"
#include "MultiGrid_Domainwall.h"
#include "MultiGrid_Domainwall-tmpl.h"

//#define USE_IMPL_IN_TMPL

// specialization on single prec.


//====================================================================

#ifdef USE_ALT_QXS
template<>
const std::string MultiGrid_Domainwall<AField_f, AField_f>::class_name
  = "MultiGrid_Domainwall< AField<float,QXS>,  AField<float,QXS> >";
#endif

#ifdef USE_ALT_ACCEL
template<>
const std::string MultiGrid_Domainwall<AField_f, AField_f>::class_name
  = "MultiGrid_Domainwall< AField<float,ACCEL>,  AField<float,ACCEL> >";
#endif


template class MultiGrid_Domainwall<AField_f, AField_f>;

//============================================================END=====
