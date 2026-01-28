/*!
      @file    fprop_Standard_lex_alt_Mixedprec.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Measurements/Fermion/fprop_Standard_lex_alt_Mixedprec.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"

#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_BiCGStab_Precond.h"


#include "lib_alt/Measurements/Fermion/fprop_Standard_lex_alt_Mixedprec-tmpl.h"

template<>
const std::string Fprop_Standard_lex_alt_Mixedprec<AField<double,SIMD2>, AField<float,SIMD2> >::class_name
 = "Fprop_Standard_lex_alt_Mixedprec<AField<double,SIMD2>, AField<float,SIMD2> >";

// class instanciation
template class Fprop_Standard_lex_alt_Mixedprec<AField<double,SIMD2>, AField<float,SIMD2> >;


//============================================================END=====
