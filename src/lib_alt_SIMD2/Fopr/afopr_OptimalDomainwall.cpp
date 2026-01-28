/*!                                                                             
        @file    afopr_OptimalDomainwall-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Fopr/afopr_OptimalDomainwall.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Tools/math_Sign_Zolotarev_Omega.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


#ifdef USE_FACTORY
namespace {
  AFopr<AField<double,SIMD2> > *create_object_with_params1(
                                             const Parameters& params)
  {
    return new AFopr_OptimalDomainwall<AField<double,SIMD2> >(params);
  }
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
                      "OptimalDomainwall", create_object_with_params1);

  AFopr<AField<float,SIMD2> > *create_object_with_params2(
                                             const Parameters& params)
  {
    return new AFopr_OptimalDomainwall<AField<float,SIMD2> >(params);
  }
  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
                      "OptimalDomainwall", create_object_with_params2);

}
#endif


#include "lib_alt/Fopr/afopr_OptimalDomainwall.h"


template class AFopr_OptimalDomainwall<AField<float,SIMD2> >;
template class AFopr_OptimalDomainwall<AField<double,SIMD2> >;

//============================================================END=====
