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

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"


// explicit instanciation.
//template class AShiftsolver_CG<AField<float,SIMD2>, AFopr<AField<float,SIMD2> > >;
template class AShiftsolver_CG<AField<double,SIMD2>, AFopr<AField<double,SIMD2> > >;

//============================================================END=====         
