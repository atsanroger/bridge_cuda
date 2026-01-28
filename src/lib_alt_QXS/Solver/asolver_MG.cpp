/*!

        @file    asolver_MG.cpp

        @brief   multigrid solver (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */
//====================================================================
#include "lib_alt/Solver/asolver_MG.h"

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

#ifdef USE_QWSLIB
#define USE_SAP_QWS
#endif

#ifdef USE_SAP_QWS
// use QWS
#include "lib_alt_QXS/Solver/asolver_SAP_QWS.h"
#endif

typedef AField<float, QXS>    AField_f;
typedef AField<double, QXS>   AField_d;

// multigrid
using MultiGrid_t = MultiGrid_Clover<AField_f, AField_f>;

// operators
using FoprD_t      = AFopr_Clover<AField_d>;
using FoprF_t      = AFopr_Clover_dd<AField_f>;
using FoprCoarse_t = AFopr_Clover_coarse<AField_f>;
#ifdef USE_SAP_QWS
#define USE_FOPR_FOR_SMOOTHER
using FoprSmoother_t = AFopr_Clover_QWS_dd<AField_f>;
#endif

// solver types
using OuterSolver_t  = ASolver_FBiCGStab<AField_d>;
using CoarseSolver_t = ASolver_BiCGStab_Cmplx<AField_f>;
#ifdef USE_SAP_QWS
using Smoother_t = ASolver_SAP_QWS<AField_f>;
#else
using Smoother_t = ASolver_SAP<AField_f>;
#define USE_SAP_FOR_SMOOTHER
#endif

#include "lib_alt/Solver/asolver_MG-tmpl.h"

template<>
const std::string ASolver_MG<AField_d>::class_name = "ASolver_MG";

template class ASolver_MG<AField_d>;
