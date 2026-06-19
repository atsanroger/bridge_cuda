/*!
      @file    afopr_Domainwall_dd.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

#include "lib_alt/Fopr/afopr_Domainwall_dd.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                      "Domainwall_dd", create_object_with_params);

  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
                      "Domainwall_dd", create_object_with_params);
}
#endif


template<>
const std::string AFopr_Domainwall_dd<AField<double,ACCEL> >::
class_name = "AFopr_Domainwall_dd<AField<double,ACCEL> >";

template<>
const std::string AFopr_Domainwall_dd<AField<float,ACCEL> >::
class_name = "AFopr_Domainwall_dd<AField<float,ACCEL> >";

#include "lib_alt/Fopr/afopr_Domainwall_dd-tmpl.h"

// class instanciation.
template class AFopr_Domainwall_dd<AField<double,ACCEL> >;
template class AFopr_Domainwall_dd<AField<float,ACCEL> >;

//============================================================END=====
