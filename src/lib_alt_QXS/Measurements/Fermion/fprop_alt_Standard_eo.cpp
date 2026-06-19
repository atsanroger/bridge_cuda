/*!
      @file    fprop_alt_Standard_eo.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
      @version $LastChangedRevision: 2562 $
*/

#define DIRECTOR_ALT_SMEARED_IMPLEMENTED

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared_alt.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/aindex_eo-inc.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string Fprop_alt_Standard_eo<AField<double, QXS> >::class_name
  = "Fprop_alt_Standard_eo<AField<double,QXS> >";


template class Fprop_alt_Standard_eo<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string Fprop_alt_Standard_eo<AField<float, QXS> >::class_name
  = "Fprop_alt_Standard_eo<AField<float,QXS> >";


template class Fprop_alt_Standard_eo<AField<float, QXS> >;

//====================================================================
// explicit instanciation for AField<half,QXS>.
#ifdef USE_QXS_FP16

template<>
const std::string Fprop_alt_Standard_eo<AField<half, QXS> >::class_name
  = "Fprop_alt_Standard_eo<AField<half,QXS> >";


template class Fprop_alt_Standard_eo<AField<half, QXS> >;

#endif
//============================================================END=====
