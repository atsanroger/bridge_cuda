/*!
      @file    fprop_alt_Standard_SAP.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-04-08 21:25:07 #$
      @version $LastChangedRevision: 2225 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_SAP.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_SAP-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.

template<>
const std::string Fprop_alt_Standard_SAP<AField<double,ACCEL> >::
class_name = "Fprop_alt_Standard_SAP<Afield<double,ACCEL> >";


template class Fprop_alt_Standard_SAP<AField<double,ACCEL> >;

//====================================================================
// explicit instanciation for AField<float,ACCEL>.

template<>
const std::string Fprop_alt_Standard_SAP<AField<float,ACCEL> >::
class_name = "Fprop_alt_Standard_SAP<Afield<float,ACCEL> >";


template class Fprop_alt_Standard_SAP<AField<float,ACCEL> >;

//====================================================================
//============================================================END=====
