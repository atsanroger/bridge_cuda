/*!
      @file    fprop_alt_Standard_eo_Richardson.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"

#include "lib_alt_Vector/Field/aindex_lex.h"
#include "lib_alt_Vector/Field/aindex_eo.h"
#include "lib_alt_Vector/Field/aindex_eo-inc.h"

#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_Richardson.h"


#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson-tmpl.h"

template<>
const std::string Fprop_alt_Standard_eo_Richardson<AField<double,VECTOR>, AField<float,VECTOR> >::class_name
  = "Fprop_alt_Standard_eo_Richardson<AField<double,VECTOR>, AField<float,VECTOR> >";

// class instanciation
template class Fprop_alt_Standard_eo_Richardson<AField<double,VECTOR>, AField<float,VECTOR> >;

//============================================================END=====
