/*!
      @file    asolver_FGMRES.cpp
      @brief   Flexible GMRES solver
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2585 $
*/
#include "lib_alt/Solver/asolver_FGMRES.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib_alt/Solver/asolver_FGMRES-tmpl.h"


// explicit instanciation.
template class ASolver_FGMRES<AField<double, QXS> >;

//============================================================END=====
