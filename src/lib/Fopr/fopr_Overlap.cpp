/*!
        @file    fopr_Overlap.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#include "Fopr/fopr_Overlap.h"
#include "Fopr/afopr_Overlap-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Overlap::register_factory();
}
#endif

template<>
const std::string AFopr_Overlap<Field>::class_name = "Fopr_Overlap";

template class AFopr_Overlap<Field>;

//============================================================END=====
