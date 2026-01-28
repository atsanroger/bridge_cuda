/*!
      @file    asolver_FBiCGStab.cpp
      @brief   Flexible BiCGStab solver
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/
#include "lib_alt/Solver/asolver_Richardson.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib_alt/Solver/asolver_Richardson-tmpl.h"


// explicit instanciation.
template class ASolver_Richardson<AField<double, QXS> >;

//============================================================END=====
