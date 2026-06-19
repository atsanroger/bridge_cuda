/*!
      @file    afopr_Staggered_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/Fopr/afopr_Staggered.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/aindex_lex.h"
#include "lib_alt_Vector/Field/afield.h"


typedef double real_t;

#include "lib_alt_Vector/inline/define_params.h"


// function template files
#include "lib_alt_Vector/Field/afield-inc.h"

// Vector functions
#include "lib_alt_Vector/BridgeVec/bridgeVec_Staggered.h"


// class template main
#include "lib_alt_Vector/Fopr/afopr_Staggered-tmpl.h"

template<>
const std::string AFopr_Staggered<AField<double,VECTOR> >::class_name
                               = "AFopr_Staggered<AField<double,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,VECTOR> >::Factory_params::Register(
                           "Staggered_eo", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Staggered<AField<double,VECTOR> >;

//============================================================END=====
