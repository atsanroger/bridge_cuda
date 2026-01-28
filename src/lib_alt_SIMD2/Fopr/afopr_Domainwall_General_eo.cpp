/*!
      @file    afopr_Domainwall_General_eo.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Fopr/afopr_Domainwall_General_eo.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
//#include "lib/Tools/math_Sign_Zolotarev_Omega.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


template<>
const std::string AFopr_Domainwall_General_eo<AField<double,SIMD2> >::
    class_name = "AFopr_Domainwall_General_eo<AField<double,SIMD2> >";

template<>
const std::string AFopr_Domainwall_General_eo<AField<float,SIMD2> >::
    class_name = "AFopr_Domainwall_General_eo<AField<float,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
                   "Domainwall_General_eo", create_object_with_params);

  bool init2 = AFopr<AField<float,SIMD2> >::Factory_params::Register(
                   "Domainwall_General_eo", create_object_with_params);
}
#endif


#include "lib_alt/Fopr/afopr_Domainwall_General_eo-tmpl.h"

// class instanciation.
template class AFopr_Domainwall_General_eo<AField<float,SIMD2> >;
template class AFopr_Domainwall_General_eo<AField<double,SIMD2> >;

//============================================================END=====
