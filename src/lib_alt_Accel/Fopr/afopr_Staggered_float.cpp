/*!
      @file    afopr_Staggered_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2535 $
*/

#ifdef USE_ACCEL_OPENACC

#include "lib_alt_Accel/Fopr/afopr_Staggered.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield.h"


typedef float real_t;

#include "lib_alt_Accel/inline/define_params.h"


// function template files
#include "lib_alt_Accel/Field/afield-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC_Staggered.h"

// class template main
#include "lib_alt_Accel/Fopr/afopr_Staggered-tmpl.h"

template<>
const std::string AFopr_Staggered<AField<float,ACCEL> >::class_name
                               = "AFopr_Staggered<AField<float,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
                           "Staggered_eo", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Staggered<AField<float,ACCEL> >;

#endif
//============================================================END=====
