/*!

        @file    asolver_MG_double.cpp

        @brief   multigrid solver (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */
//====================================================================
#include "lib_alt/Solver/asolver_MG_double.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h" // matching btw coasre/fine lattices

// clover specific
#include "lib_alt_QXS/Fopr/afopr_Clover.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_coarse.h"
#include "lib_alt/Solver/MultiGrid_Clover.h"


#include "lib_alt/Solver/asolver_SAP_MINRES.h"
#include "lib_alt/Solver/asolver_SAP.h"


//typedef AField<float, QXS> AField_f;
typedef AField<double, QXS> AField_d;

// multigrid
using MultiGrid_t = MultiGrid_Clover<AField_d, AField_d>;

// operators
using FoprD_t      = AFopr_Clover<AField_d>;
using FoprF_t      = AFopr_Clover_dd<AField_d>;
using FoprCoarse_t = AFopr_Clover_coarse<AField_d>;

// solver types
using OuterSolver_t  = ASolver_FBiCGStab<AField_d>;
using CoarseSolver_t = ASolver_BiCGStab_Cmplx<AField_d>;
//using CoarseSolver_t = ASolver_BiCGStab< AField_f>;
using Smoother_t = ASolver_SAP<AField_d>;
#define USE_SAP_FOR_SMOOTHER


#include "lib_alt/Solver/asolver_MG_double-tmpl.h"

template<>
const std::string ASolver_MG_double<AField_d>::class_name = "ASolver_MG_double";

template class ASolver_MG_double<AField_d>;
