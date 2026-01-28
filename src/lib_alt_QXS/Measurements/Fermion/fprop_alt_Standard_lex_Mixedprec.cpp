/*!
      @file    fprop_alt_Standard_lex_Mixedprec.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_FBiCGStab.h"
#include "lib_alt/Solver/asolver_BiCGStab_Precond.h"


#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec-tmpl.h"

template<>
const std::string Fprop_alt_Standard_lex_Mixedprec<AField<double, QXS>, AField<float, QXS> >::class_name
  = "Fprop_alt_Standard_lex_Mixedprec<AField<double,QXS>, AField<float,QXS> >";

// class instanciation
template class Fprop_alt_Standard_lex_Mixedprec<AField<double, QXS>, AField<float, QXS> >;


//============================================================END=====
