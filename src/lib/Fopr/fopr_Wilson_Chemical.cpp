/*!
        @file    fopr_Wilson_Chemical.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/fopr_Wilson_Chemical.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Wilson_Chemical::register_factory();
}
#endif


#include "Fopr/afopr_Wilson_Chemical-tmpl.h"

template<>
const std::string Fopr_Wilson_Chemical::class_name
  = "Fopr_Wilson_Chemical";

template class AFopr_Wilson_Chemical<Field>;

//============================================================END=====
