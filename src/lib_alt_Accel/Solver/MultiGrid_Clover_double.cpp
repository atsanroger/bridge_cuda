/*!
        @file    MultiGrid_Clover_double.cpp
        @brief   MultiGrid operation for Clover fermion (Accel version)
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2229 $
 */

//====================================================================

#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Fopr/afopr_Clover.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"

// template for MultiGrid_Clover
#include "lib_alt/Solver/MultiGrid_Clover.h"
#include "lib_alt/Solver/MultiGrid_Clover-tmpl.h"

typedef AField<double, ACCEL> AField_d;

template<>
const std::string MultiGrid_Clover< AField_d, AField_d >::class_name
 = "MultiGrid_Clover< AField<double,ACCEL>,  AField<double,ACCEL> >";
template class MultiGrid_Clover< AField_d, AField_d >;

//============================================================END=====
