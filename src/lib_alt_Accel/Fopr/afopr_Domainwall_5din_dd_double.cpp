/*!
      @file    afopr_Domainwall_5din_dd_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"

// C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

// Bridge++ core library header files
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"


typedef double real_t;

#include "lib_alt_Accel/inline/define_params.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Wilson.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Domainwall.h"


// template file
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd-tmpl.h"

template<>
const std::string AFopr_Domainwall_5din_dd<AField<double, ACCEL> >
::class_name = "AFopr_Domainwall_5din_dd<AField<double,ACCEL> >";


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, ACCEL> >::Factory_params::Register(
    "Domainwall_5din_dd", create_object_with_params1);
}
#endif


// explicit instanciation
template class AFopr_Domainwall_5din_dd<AField<double, ACCEL> >;

//============================================================END=====
