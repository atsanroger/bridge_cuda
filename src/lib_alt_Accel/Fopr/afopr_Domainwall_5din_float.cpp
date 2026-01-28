/*!
      @file    afopr_Domainwall_5din_float.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2553 $
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"

// C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

// Bridge++ core library header files
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"


typedef float real_t;

#include "lib_alt_Accel/inline/define_params.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Domainwall.h"


// template file
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din-tmpl.h"

template<>
const std::string AFopr_Domainwall_5din<AField<float, ACCEL> >
::class_name = "AFopr_Domainwall_5din<AField<float,ACCEL> >";


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<float, ACCEL> >::Factory_params::Register(
    "Domainwall_5din", create_object_with_params1);
}
#endif


// explicit instanciation
template class AFopr_Domainwall_5din<AField<float, ACCEL> >;

//============================================================END=====
