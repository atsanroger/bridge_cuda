/*!
      @file    asolver_MG.cpp
      @brief   multigrid solver (Accel version)
      @author  KANAMORI Issaku (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate::  #$
      @version $LastChangedRevision: 2227 $
 */

//====================================================================
#include "lib_alt/Solver/asolver_MG.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"

// clover specific
#include "lib_alt_Accel/Fopr/afopr_Clover.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_coarse.h"
#include "lib_alt/Solver/MultiGrid_Clover.h"


#include "lib_alt/Solver/asolver_SAP_MINRES.h"
#include "lib_alt/Solver/asolver_SAP.h"


typedef AField<float, ACCEL> AField_f;
typedef AField<double, ACCEL> AField_d;

// multigrid
using MultiGrid_t = MultiGrid_Clover<AField_f, AField_f>;

// operators
using FoprD_t = AFopr_Clover<AField_d>;
using FoprF_t = AFopr_Clover_dd<AField_f>;
using FoprCoarse_t = AFopr_Clover_coarse<AField_f >;

// solver types
using OuterSolver_t  = ASolver_FBiCGStab<AField_d>;
using CoarseSolver_t = ASolver_BiCGStab_Cmplx< AField_f>;
//using CoarseSolver_t = ASolver_BiCGStab< AField_f>;
using Smoother_t     = ASolver_SAP<AField_f>;
#define USE_SAP_FOR_SMOOTHER


#include "lib_alt/Solver/asolver_MG-tmpl.h"

template<>
const std::string ASolver_MG<AField_d>::class_name = "ASolver_MG";

template class ASolver_MG<AField_d>;

