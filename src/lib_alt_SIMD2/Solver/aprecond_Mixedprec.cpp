/*!
      @file    aprecond_Mixedprec.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include <cassert>

#include "lib_alt/Solver/aprecond_Mixedprec.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/index_eo_alt.h"


#include "lib_alt/Solver/aprecond_Mixedprec-tmpl.h"


// explicit instanciation.
template class APrecond_Mixedprec<AField<double,SIMD2>, AField<float,SIMD2> >;

//============================================================END=====
