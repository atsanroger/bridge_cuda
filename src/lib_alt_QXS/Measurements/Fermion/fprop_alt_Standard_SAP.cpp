/*!
      @file    fprop_alt_Standard_SAP.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_SAP.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_SAP-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.

template<>
const std::string Fprop_alt_Standard_SAP<AField<double, QXS> >::
class_name = "Fprop_alt_Standard_SAP<Afield<double,QXS> >";


template class Fprop_alt_Standard_SAP<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.

template<>
const std::string Fprop_alt_Standard_SAP<AField<float, QXS> >::
class_name = "Fprop_alt_Standard_SAP<Afield<float,QXS> >";


template class Fprop_alt_Standard_SAP<AField<float, QXS> >;

//============================================================END=====
