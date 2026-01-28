/*!
        @file    fopr_Smeared_alt.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-11-13 13:02:59 #$

        @version $LastChangedRevision: 2553 $
*/

#ifdef USE_ALT_CODE

#include "Fopr/fopr_Smeared_alt.h"
#include "Fopr/afopr_Smeared_alt-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Smeared_alt::register_factory();
}
#endif
template<>
const std::string AFopr_Smeared_alt<Field>::class_name = "Fopr_Smeared_alt";

template class AFopr_Smeared_alt<Field>;

#endif

//============================================================END=====
