/*!

        @file    $Id: MultiGrid_Clover.cpp #$

        @brief   MultiGrid operation for Clover fermion (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */

//====================================================================

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
#include "lib_alt_QXS/Fopr/afopr_Clover.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"

// template for MultiGrid_Clover
#include "lib_alt/Solver/MultiGrid_Clover.h"
#include "lib_alt/Solver/MultiGrid_Clover-tmpl.h"

typedef AField<float, QXS> AField_f;
//#define USE_IMPL_IN_TMPL

// specialization on single prec.


//====================================================================
template<>
const std::string MultiGrid_Clover<AField_f, AField_f>::class_name = "MultiGrid_Clover< AField<float,QXS>,  AField<float,QXS> >";
template class MultiGrid_Clover<AField_f, AField_f>;


//============================================================END=====
