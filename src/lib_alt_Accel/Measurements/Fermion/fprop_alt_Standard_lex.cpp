/*!
      @file    fprop_alt_Standard_lex.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string Fprop_alt_Standard_lex<AField<double,ACCEL> >::class_name
                         = "Fprop_alt_Standard_lex<Afield<double,ACCEL> >";


template class Fprop_alt_Standard_lex<AField<double,ACCEL> >;

//====================================================================
// explicit instanciation for AField<float,ACCEL>.
template<>
const std::string Fprop_alt_Standard_lex<AField<float,ACCEL> >::class_name
                          = "Fprop_alt_Standard_lex<Afield<float,ACCEL> >";


template class Fprop_alt_Standard_lex<AField<float,ACCEL> >;

//============================================================END=====
