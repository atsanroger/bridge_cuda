/*!
      @file    afopr_Wilson_eo_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2535 $
*/

#ifdef USE_ACCEL_OPENACC

#include "lib_alt_Accel/Fopr/afopr_Wilson_eo.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/aindex_eo.h"


typedef double real_t;

#include "lib_alt_Accel/inline/define_params.h"

// function template files
#include "lib_alt_Accel/Field/afield-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Wilson.h"

// class template
#include "lib_alt_Accel/Fopr/afopr_Wilson_eo-tmpl.h"

template<>
const std::string AFopr_Wilson_eo<AField<double,ACCEL> >::class_name
                             = "AFopr_Wilson_eo<AField<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                           "Wilson_eo", create_object_with_params);
}
#endif


// explicit instanciation.
template class AFopr_Wilson_eo<AField<double,ACCEL> >;

#endif
//============================================================END=====
