/*!
      @file    asolver_BiCGStab_Precond.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/


#include "lib_alt/Solver/asolver_BiCGStab_Precond.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_BiCGStab_Precond-tmpl.h"


// explicit instanciation.
template class ASolver_BiCGStab_Precond<AField<double,ACCEL> >;

//============================================================END=====
