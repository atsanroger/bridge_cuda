/*!
      @file    afopr_Domainwall_5din_eo_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

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
#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Domainwall_SU3_double-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_qxs_parts-inc.h"


// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Domainwall.h"

// template file
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo-tmpl.h"


template<>
const std::string AFopr_Domainwall_5din_eo<AField<double, QXS> >
::class_name = "AFopr_Domainwall_5din_eo<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Domainwall_5din_eo", create_object_with_params1);
  init1 &= AFopr<AField<double, QXS> >::Factory_params::Register(
    "Domainwall_General_5din_eo", create_object_with_params1);
}
#endif


// explicit instanciation
template class AFopr_Domainwall_5din_eo<AField<double, QXS> >;



//============================================================END=====
