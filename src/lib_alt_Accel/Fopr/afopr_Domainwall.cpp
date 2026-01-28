/*!
      @file    afopr_Domainwall.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib/Fopr/afopr_Domainwall.h"

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
                      "Domainwall", create_object_with_params);

  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
                      "Domainwall", create_object_with_params);
}
#endif


template<>
const std::string AFopr_Domainwall<AField<double,ACCEL> >::
class_name = "AFopr_Domainwall<AField<double,ACCEL> >";

template<>
const std::string AFopr_Domainwall<AField<float,ACCEL> >::
class_name = "AFopr_Domainwall<AField<float,ACCEL> >";

#include "lib/Fopr/afopr_Domainwall-tmpl.h"

// class instanciation.
template class AFopr_Domainwall<AField<double,ACCEL> >;
template class AFopr_Domainwall<AField<float,ACCEL> >;

//============================================================END=====
