/*!
      @file    afopr_Domainwall_PVprec_float.cpp
      @brief   float instantiation + factory registration of
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


typedef float real_t;

#include "lib_alt_Accel/inline/define_params.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"


template<>
const std::string AFopr_Domainwall_PVprec<AField<float, ACCEL> >
::class_name = "AFopr_Domainwall_PVprec<AField<float,ACCEL> >";


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr_Domainwall_PVprec<AField<float, ACCEL> >::register_factory();
}
#endif


// explicit instanciation
template class AFopr_Domainwall_PVprec<AField<float, ACCEL> >;

//============================================================END=====
