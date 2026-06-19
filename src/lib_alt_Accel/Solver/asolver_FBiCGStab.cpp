/*!
      @file    asolver_FBiCGStab.cpp
      @brief   Flexible BiCGStab solver
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2192 $
*/
#include "lib_alt/Solver/asolver_FBiCGStab.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt/Solver/asolver_FBiCGStab-tmpl.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

// explicit instanciation.
template class ASolver_FBiCGStab<AField<double,ACCEL> >;

//============================================================END=====
