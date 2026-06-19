/*!
      @file    afopr_Domainwall_5din_eo_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2613 $
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo.h"

// C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

// Bridge++ core library header files
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_eo.h"

typedef double real_t;

#include "lib_alt_Accel/inline/define_params.h"

#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_eo-inc.h"

#include "lib_alt_Accel/BridgeACC/bridgeACC.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Domainwall.h"


// template file
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, ACCEL> >::Factory_params::Register(
    "Domainwall_5din_eo", create_object_with_params1);
}
#endif


// explicit instanciation
template<>
const std::string AFopr_Domainwall_5din_eo<AField<double, ACCEL> >
::class_name = "AFopr_Domainwall_5din_eo<AField<double,ACCEL> >";

template class AFopr_Domainwall_5din_eo<AField<double, ACCEL> >;

//============================================================END=====
