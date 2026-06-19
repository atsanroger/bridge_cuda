/*!

        @file    asolver_MG_dw.cpp

        @brief   multigrid solver for domainwall fermion (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2586 $

 */
//====================================================================
#include "lib_alt/Solver/asolver_MG.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Tools/randomNumberManager.h"


#ifdef USE_ALT_QXS
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"

// Domainwall specific
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_coarse.h"

typedef AField<float, QXS>    AField_f;
typedef AField<double, QXS>   AField_d;
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"

// Domainwall specific
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_coarse.h"

typedef AField<float,  ACCEL>  AField_f;
typedef AField<double, ACCEL>  AField_d;
#endif

#include "MultiGrid_Domainwall.h"

#include "lib_alt/Solver/asolver_Richardson.h"
#include "lib_alt/Solver/asolver_GMRES_m_Cmplx.h"
#include "lib_alt/Solver/asolver_FGMRES.h"
#include "lib_alt/Solver/asolver_CGNR.h"
#include "asolver_SAP_MINRES_dw.h"
#include "asolver_SAP_MINRES_dw2.h"
#include "asolver_SAP_dw.h"
#include "asolver_SAP_dw2.h"

// multigrid
using MultiGrid_t = MultiGrid_Domainwall<AField_f, AField_f>;

// operators
using FoprD_t      = AFopr_Domainwall_5din<AField_d>;
using FoprF_t      = AFopr_Domainwall_5din_dd<AField_f>;
using FoprCoarse_t = AFopr_Domainwall_coarse<AField_f>;
// solver types
//using OuterSolver_t  = ASolver_FBiCGStab<AField_d>;
using OuterSolver_t  = ASolver_FGMRES<AField_d>;
//using OuterSolver_t  = ASolver_Richardson<AField_d>;

using CoarseSolver_t = ASolver_BiCGStab_Cmplx<AField_f>;
//using Smoother_t = ASolver_SAP_dw<AField_f>;
//using Smoother_t = ASolver_SAP_dw2<AField_f>;
//#define USE_SAP_FOR_SMOOTHER
////using Smoother_t = ASolver_GMRES_m_Cmplx<AField_f>;
using Smoother_t = ASolver_CGNR<AField_f>;


#include "asolver_MG_dw-tmpl.h"

template<>
const std::string ASolver_MG_dw<AField_d>::class_name = "ASolver_MG_dw<AField<double, QXS> >";

template class ASolver_MG_dw<AField_d>;
