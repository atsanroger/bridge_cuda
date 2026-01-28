/*!
        @file    fopr_Domainwall.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#include "Fopr/fopr_Domainwall.h"
#include "Fopr/afopr_Domainwall-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Domainwall::register_factory();
}
#endif
template<>
const std::string AFopr_Domainwall<Field>::class_name = "Fopr_Domainwall";

template class AFopr_Domainwall<Field>;


//============================================================END=====
