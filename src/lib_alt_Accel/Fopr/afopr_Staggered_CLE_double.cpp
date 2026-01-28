/*!
      @file    afopr_Staggered_CLE_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2535 $
*/

#ifdef USE_ACCEL_OPENACC


#include "lib_alt_Accel/Fopr/afopr_Staggered_CLE.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield.h"

typedef double real_t;

// function template files
#include "lib_alt_Accel/Field/afield-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC_Staggered.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField_Gauge.h"

// class template main
#include "lib_alt_Accel/Fopr/afopr_Staggered_CLE-tmpl.h"

template<>
const std::string AFopr_Staggered_CLE<AField<double,ACCEL> >::class_name
                      = "AFopr_Staggered_CLE<AField<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                           "Staggered_CLE", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Staggered_CLE<AField<double,ACCEL> >;


#endif
//============================================================END=====
