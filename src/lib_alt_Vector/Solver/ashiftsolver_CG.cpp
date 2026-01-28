/*!
      @file    ashiftsolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib/Solver/ashiftsolver_CG.h"
#include "lib/Solver/ashiftsolver_CG-tmpl.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt_Vector/Field/afield.h"
#include "lib_alt_Vector/Field/afield-inc.h"


// explicit instanciation.
//template class AShiftsolver_CG<AField<float,VECTOR>, AFopr<AField<float,VECTOR> > >;
template class AShiftsolver_CG<AField<double,VECTOR>, AFopr<AField<double,VECTOR> > >;

//============================================================END=====         
