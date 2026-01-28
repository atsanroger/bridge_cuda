/*!
      @file    fprop_alt_Standard_lex.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"
#include "lib_alt_Vector/Field/aindex_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,VECTOR>.
template<>
const std::string Fprop_alt_Standard_lex<AField<double,VECTOR> >::class_name
                         = "Fprop_alt_Standard_lex<Afield<double,VECTOR> >";


template class Fprop_alt_Standard_lex<AField<double,VECTOR> >;

//====================================================================
// explicit instanciation for AField<float,VECTOR>.
template<>
const std::string Fprop_alt_Standard_lex<AField<float,VECTOR> >::class_name
                          = "Fprop_alt_Standard_lex<Afield<float,VECTOR> >";


template class Fprop_alt_Standard_lex<AField<float,VECTOR> >;

//============================================================END=====
