/*!
      @file    ashiftsolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib/Solver/ashiftsolver_CG.h"
#include "lib/Solver/ashiftsolver_CG-tmpl.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"


// explicit instanciation.
//template class AShiftsolver_CG<AField<float,ACCEL>, AFopr<AField<float,ACCEL> > >;
template class AShiftsolver_CG<AField<double,ACCEL>, AFopr<AField<double,ACCEL> > >;

//============================================================END=====         
