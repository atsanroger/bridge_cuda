/*!
      @file    afopr_Domainwall_5din_eo_half.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/
#ifdef USE_FP16

#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo.h"

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
#define  VLEN     VLENH
#define  VLENX    VLENXH
#define  VLENY    VLENYH

typedef half real_t;

#include "lib_alt_QXS/inline/vsimd_half-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_half-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_half-inc.h"
#include "lib_alt_QXS/inline/vsimd_Domainwall_SU3_half-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_qxs_parts-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Domainwall.h"

// template file
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo-tmpl.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  AFopr<AField<half, QXS> > *create_object_with_params1(
    const Parameters& params)
  {
    return new AFopr_Domainwall_5din_eo<AField<half, QXS> >(params);
  }


  bool init1 = AFopr<AField<half, QXS> >::Factory_params::Register(
    "Domainwall_5din_eo", create_object_with_params1);
  // temporal for transition of the name
  init1 &= AFopr<AField<half, QXS> >::Factory_params::Register(
    "Domainwall_General_5din_eo", create_object_with_params1);
}
#endif

template<>
const std::string AFopr_Domainwall_5din_eo<AField<half, QXS> >
::class_name = "AFopr_Domainwall_5din_eo<AField<half,QXS> >";

// explicit instanciation
template class AFopr_Domainwall_5din_eo<AField<half, QXS> >;

#endif
//============================================================END=====
