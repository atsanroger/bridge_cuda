/*!
      @file    fprop_alt_Standard_lex_Mixedprec.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"
#include "lib_alt_Vector/Field/aindex_lex.h"

#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_FBiCGStab.h"
#include "lib_alt/Solver/asolver_BiCGStab_Precond.h"


#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec-tmpl.h"

template<>
const std::string Fprop_alt_Standard_lex_Mixedprec<AField<double,VECTOR>, AField<float,VECTOR> >::class_name
 = "Fprop_alt_Standard_lex_Mixedprec<AField<double,VECTOR>, AField<float,VECTOR> >";

// class instanciation
template class Fprop_alt_Standard_lex_Mixedprec<AField<double,VECTOR>, AField<float,VECTOR> >;


//============================================================END=====
