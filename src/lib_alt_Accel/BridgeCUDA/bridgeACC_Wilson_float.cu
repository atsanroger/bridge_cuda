/*!
      @file    bridgeACC_Wilson_float.cu
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


namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "src/mult_Wilson_cuda-inc.h"
#include "src/mult_Wilson_cuda2-inc.h"
#include "src/mult_Wilson_dd_cuda-inc.h"

// QDW Stubs for float
void mult_wilson_qdw_D_dirac(float* v2, float* u, float* v1,
                             int* Nsize, int* bc, float kappa)
{
  printf("Error: QDW mult only supported for double precision.\n");
  exit(1);
}

void mult_wilson_qdw_D_chiral(float* v2, float* u, float* v1,
                              int* Nsize, int* bc, float kappa)
{
  printf("Error: QDW mult only supported for double precision.\n");
  exit(1);
}


}

//============================================================END=====
