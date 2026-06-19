/*!
      @file    afopr_Domainwall.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib/Fopr/afopr_Domainwall.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Tools/math_Sign_Zolotarev_Omega.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"

#define IMPL VECTOR

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,IMPL> >::Factory_params::Register(
                      "Domainwall", create_object_with_params);

  bool init2 = AFopr<AField<float,IMPL> >::Factory_params::Register(
                      "Domainwall", create_object_with_params);
}
#endif


template<>
const std::string AFopr_Domainwall<AField<double,IMPL> >::class_name
                       = "AFopr_Domainwall<AField<double,IMPL> >";

template<>
const std::string AFopr_Domainwall<AField<float,IMPL> >::class_name
                       = "AFopr_Domainwall<AField<float,IMPL> >";

#include "lib/Fopr/afopr_Domainwall-tmpl.h"

// class instanciation.
template class AFopr_Domainwall<AField<float,IMPL> >;
template class AFopr_Domainwall<AField<double,IMPL> >;

//============================================================END=====
