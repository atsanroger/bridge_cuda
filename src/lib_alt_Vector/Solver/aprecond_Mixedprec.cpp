/*!
      @file    aprecond_Mixedprec.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include <cassert>

#include "lib_alt/Solver/aprecond_Mixedprec.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"
#include "lib_alt_Vector/Field/aindex_lex.h"
#include "lib_alt_Vector/Field/aindex_eo.h"


#include "lib_alt/Solver/aprecond_Mixedprec-tmpl.h"


// explicit instanciation.
template class APrecond_Mixedprec<AField<double,VECTOR>, AField<float,VECTOR> >;

//============================================================END=====
