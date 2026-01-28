/*!
      @file    asolver_FGMRES.cpp
      @brief   Flexible GMRES solver
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2585 $
*/
#include "lib_alt/Solver/asolver_FGMRES.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib_alt/Solver/asolver_FGMRES-tmpl.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"


// explicit instanciation.
template class ASolver_FGMRES<AField<double, ACCEL> >;

//============================================================END=====
