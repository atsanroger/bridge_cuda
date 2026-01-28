/*!
      @file    fprop_alt_Standard_lex.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
      @version $LastChangedRevision: 2562 $
*/

#define DIRECTOR_ALT_SMEARED_IMPLEMENTED

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared_alt.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string Fprop_alt_Standard_lex<AField<double, QXS> >::class_name
  = "Fprop_alt_Standard_lex<Afield<double,QXS> >";


template class Fprop_alt_Standard_lex<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string Fprop_alt_Standard_lex<AField<float, QXS> >::class_name
  = "Fprop_alt_Standard_lex<Afield<float,QXS> >";


template class Fprop_alt_Standard_lex<AField<float, QXS> >;

//====================================================================
// explicit instanciation for AField<half,QXS>.
#ifdef USE_QXS_FP16

template<>
const std::string Fprop_alt_Standard_lex<AField<half, QXS> >::class_name
  = "Fprop_alt_Standard_lex<Afield<half,QXS> >";


template class Fprop_alt_Standard_lex<AField<half, QXS> >;

#endif
//============================================================END=====
