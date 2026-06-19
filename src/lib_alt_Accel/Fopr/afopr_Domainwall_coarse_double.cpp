/*!
      @file    afopr_Domainwall_coarse_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2223 $
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_coarse.h"

#include "lib/ResourceManager/threadManager.h"
#include "complexTraits.h"

#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;
typedef typename ComplexTraits<double>::complex_t complex_t;


#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"

#include "lib_alt_Accel/Fopr/afopr_common_th-inc.h"


// template definition
#include "lib_alt_Accel/Fopr/afopr_Domainwall_coarse-tmpl.h"

template<>
const std::string AFopr_Domainwall_coarse<AField<double,ACCEL> >::class_name
                      = "AFopr_Domainwall_coarse<AField<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                          "Domainwall_coarse", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Domainwall_coarse<AField<double,ACCEL> >;

//============================================================END=====
