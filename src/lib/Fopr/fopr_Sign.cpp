/*!
        @file    fopr_Sign.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#include "Fopr/fopr_Sign.h"
#include "Fopr/afopr_Sign-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Sign::register_factory();
}
#endif

template<>
const std::string AFopr_Sign<Field>::class_name = "Fopr_Sign";

template class AFopr_Sign<Field>;


//============================================================END=====
