/*!
      @file    afopr_Domainwall_5din_dd_double.cpp
      @brief
      @author  Issaku Kanamori
               $LastChangedBy: $
      @date    $LastChangedDate: $
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd.h"

// C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

// Bridge++ core library header files
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

// vector length
#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_Domainwall_SU3_float-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
//#include "lib_alt_QXS/Fopr/mult_Wilson_qxs_parts-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Domainwall.h"
#include "lib_alt_QXS/BridgeQXS/src/mult_common_parts_qxs-inc.h"
#include "lib_alt_QXS/BridgeQXS/src/mult_Wilson_parts_qxs-inc.h"

// template file
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  AFopr<AField<float, QXS> > *create_object_with_params1(
    const Parameters& params)
  {
    return new AFopr_Domainwall_5din_dd<AField<float, QXS> >(params);
  }


  bool init1 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Domainwall_5din_dd", create_object_with_params1);
}
#endif

// explicit instanciation
template<>
const std::string AFopr_Domainwall_5din_dd<AField<float, QXS> >
::class_name = "AFopr_Domainwall_5din_dd<AField<float,QXS> >";

template class AFopr_Domainwall_5din_dd<AField<float, QXS> >;

//============================================================END=====
