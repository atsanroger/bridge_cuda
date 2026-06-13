/*!

        @file    asolver_MG_dw.cpp

        @brief   multigrid solver for domainwall fermion (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2595 $

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
#include "lib_alt_Accel/Fopr/afopr_Domainwall_PVprec.h"

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
#include "asolver_PVexact_dw.h"

#include <cstdio>     // fflush
#include <unistd.h>   // dup, dup2, close (mute fd 1 around the setup smoother)
#include <fcntl.h>    // open, O_WRONLY

// multigrid
using MultiGrid_t = MultiGrid_Domainwall<AField_f, AField_f>;

// operators
// Option 2 (A-multigrid): the OUTER double operator is the positive, well-
// conditioned target A = (D_PV C_PV^{-1})^dag C^{-1} D (the paper's MG target),
// NOT the bare non-normal D (on which GMRES/FGMRES diverge). The outer solve
// becomes exactly the GMRES-on-A we validated in the spectrum test.
using FoprD_t      = AFopr_Domainwall_PVprec<AField_d>;
using FoprF_t      = AFopr_Domainwall_5din_dd<AField_f>;
using FoprCoarse_t = AFopr_Domainwall_coarse<AField_f>;
// solver types
// Option 2, step 2 (A-preconditioner): the outer is FGMRES on the positive A,
// preconditioned by the clean two-grid V-cycle on A (smoother=GMRES-on-A +
// A-Galerkin coarse). FGMRES (flexible) is required because the V-cycle is a
// non-constant preconditioner. The consistency baseline (no precond, ~2221
// double iters) is the GMRES_m_Cmplx line kept below for comparison.
//using OuterSolver_t  = ASolver_FBiCGStab<AField_d>;
using OuterSolver_t  = ASolver_FGMRES<AField_d>;
//using OuterSolver_t  = ASolver_GMRES_m_Cmplx<AField_d>;
//using OuterSolver_t  = ASolver_Richardson<AField_d>;

//using CoarseSolver_t = ASolver_BiCGStab_Cmplx<AField_f>;
using CoarseSolver_t = ASolver_GMRES_m_Cmplx<AField_f>;  // V2: GMRES coarse (BiCGStab breaks down -> nan)
// SAP smoother (DD-alphaAMG for DWF): block SAP solve + PV solver.
// CGNR-on-bare-D smoother diverges at small mq/M0=1.0 (near-critical D^dag D),
// which is exactly the DWF pathology the PV-preconditioned MG is meant to avoid.
// Sweep binaries prebuilt in sweep_bins/: v1=SAP+BiCGStab, v2=SAP+GMRES,
// v3=GMRES smoother+GMRES coarse, v4=CGNR smoother+GMRES coarse.
// Source left in the recommended state: SAP smoother + GMRES coarse (=v2).
#define USE_SAP_FOR_SMOOTHER
using Smoother_t = ASolver_SAP_dw<AField_f>;
//using Smoother_t = ASolver_SAP_dw2<AField_f>;
//using Smoother_t = ASolver_GMRES_m_Cmplx<AField_f>;
//using Smoother_t = ASolver_CGNR<AField_f>;


#include "asolver_MG_dw-tmpl.h"

template<>
#ifdef USE_ALT_ACCEL
const std::string ASolver_MG_dw<AField_d>::class_name = "ASolver_MG_dw<AField<double, ACCEL> >";
#else
const std::string ASolver_MG_dw<AField_d>::class_name = "ASolver_MG_dw<AField<double, QXS> >";
#endif

template class ASolver_MG_dw<AField_d>;
