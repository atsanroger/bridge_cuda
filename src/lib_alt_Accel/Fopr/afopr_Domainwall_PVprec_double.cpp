/*!
      @file    afopr_Domainwall_PVprec_double.cpp
      @brief   double instantiation + factory registration of
               AFopr_Domainwall_PVprec.
      @author  Wei-Lun Chen
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_PVprec.h"

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


template<>
const std::string AFopr_Domainwall_PVprec<AField<double, ACCEL> >
::class_name = "AFopr_Domainwall_PVprec<AField<double,ACCEL> >";


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr_Domainwall_PVprec<AField<double, ACCEL> >::register_factory();
}
#endif


// explicit instanciation
template class AFopr_Domainwall_PVprec<AField<double, ACCEL> >;

//============================================================END=====
