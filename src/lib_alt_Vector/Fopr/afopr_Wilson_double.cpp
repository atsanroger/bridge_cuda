/*!
      @file    afopr_Wilson_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/Fopr/afopr_Wilson.h"

#include "lib/ResourceManager/threadManager.h"

typedef double real_t;

#include "lib_alt_Vector/inline/define_params.h"
//#include "lib_alt_Vector/inline/mult_th-inc.h"

#include "lib_alt_Vector/Field/aindex_lex.h"
#include "lib_alt_Vector/Field/afield.h"

// function template files
#include "lib_alt_Vector/Field/afield-inc.h"

// Vector functions header files
#include "lib_alt_Vector/BridgeVec/bridgeVec_afield.h"
#include "lib_alt_Vector/BridgeVec/bridgeVec_Wilson.h"

// class template main
#include "lib_alt_Vector/Fopr/afopr_Wilson-tmpl.h"

template<>
const std::string AFopr_Wilson<AField<double,VECTOR> >::class_name
                            = "AFopr_Wilson<AField<double,VECTOR> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,VECTOR> >::Factory_params::Register(
                               "Wilson", create_object_with_params);
}
#endif


// explicit instanciation
template class AFopr_Wilson<AField<double,VECTOR> >;

//============================================================END=====
