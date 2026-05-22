/*!
      @file    bridgeACC_Domainwall_5din_float.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"
#include "bridgeACC.h"
#include "bridgeACC_AField.h"


typedef float real_t;
#include "inline/mult_Wilson_cuda_inline-inc.h"
//#include "src/mult_Wilson_cuda2-inc.h"


namespace BridgeACC {

#define SU3_3RD_ROW_RECONST
//#define USE_5D_HOPB_V2  // shared-memory gauge link sharing for 5D EO hopb kernel
//#define USE_LUINV_V2    // thread-local vt_fwd[NS][ND] for LU(dag)inv (drops vp round-trip)

#include "src/mult_Domainwall_5din_cuda-inc.h"
#include "src/mult_Domainwall_5din_eo_cuda-inc.h"
#include "src/mult_Domainwall_5din_dd_cuda-inc.h"
#include "src/mult_Domainwall_5din_cuda2-inc.h"
#include "src/mult_Domainwall_5din_eo_inv_cuda-inc.h"

// QDW EO stubs — QDW arithmetic is double-only; print error if called for float
void mult_domainwall_5din_ee_5dir_dirac_qdw(
    float*, float*, float, float, int, float*, float*, float, int*)
{ printf("Error: QDW mult only supported for double precision.\n"); exit(1); }

void mult_domainwall_5din_eo_5dir_dirac_qdw(
    float*, float*, float, float, int, float*, float*, float, int*)
{ printf("Error: QDW mult only supported for double precision.\n"); exit(1); }

void mult_domainwall_5din_ee_5dirdag_dirac_qdw(
    float*, float*, float, float, int, float*, float*, float, int*)
{ printf("Error: QDW mult only supported for double precision.\n"); exit(1); }

void mult_domainwall_5din_eo_5dirdag_dirac_qdw(
    float*, float*, float, float, int, float*, float*, float, int*)
{ printf("Error: QDW mult only supported for double precision.\n"); exit(1); }

void mult_domainwall_5din_eo_hopb_qdw_dirac_5d(
    float*, float*, float*, int, int*, int*, int*, int, int, int)
{ printf("Error: QDW mult only supported for double precision.\n"); exit(1); }

void mult_domainwall_5din_ee_LUinv_dirac_qdw(
    float*, float*, int, int*, float)
{ printf("Error: QDW LU inv only supported for double precision.\n"); exit(1); }

void mult_domainwall_5din_ee_LUdaginv_dirac_qdw(
    float*, float*, int, int*, float)
{ printf("Error: QDW LUdag inv only supported for double precision.\n"); exit(1); }

}

//============================================================END=====
