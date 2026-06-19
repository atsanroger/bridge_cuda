/*!
      @file    afopr_Domainwall.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib/Fopr/afopr_Domainwall.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Domainwall", create_object_with_params);

  bool init2 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Domainwall", create_object_with_params);
}
#endif


template<>
const std::string AFopr_Domainwall<AField<double, QXS> >::class_name
  = "AFopr_Domainwall<AField<double,QXS> >";

template<>
const std::string AFopr_Domainwall<AField<float, QXS> >::class_name
  = "AFopr_Domainwall<AField<float,QXS> >";


#include "lib/Fopr/afopr_Domainwall-tmpl.h"

// class instanciation.
template class AFopr_Domainwall<AField<double, QXS> >;
template class AFopr_Domainwall<AField<float, QXS> >;

//============================================================END=====
