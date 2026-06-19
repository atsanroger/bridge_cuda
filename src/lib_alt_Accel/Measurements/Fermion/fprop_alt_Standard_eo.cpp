/*!
      @file    fprop_alt_Standard_eo.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-29 15:46:01 #$
      @version $LastChangedRevision: 2542 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/aindex_eo.h"
#include "lib_alt_Accel/Field/aindex_eo-inc.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string Fprop_alt_Standard_eo<AField<double,ACCEL> >::class_name
                          = "Fprop_alt_Standard_eo<Afield<double,ACCEL> >";


template class Fprop_alt_Standard_eo<AField<double,ACCEL> >;

//====================================================================
// explicit instanciation for AField<float,ACCEL>.
template<>
const std::string Fprop_alt_Standard_eo<AField<float,ACCEL> >::class_name
                           = "Fprop_alt_Standard_eo<Afield<float,ACCEL> >";


template class Fprop_alt_Standard_eo<AField<float,ACCEL> >;

//============================================================END=====
