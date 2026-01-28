/*!
      @file    aprecond.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Solver/aprecond.h"

#include "lib_alt_Vector/Field/afield.h"

// explicit instanciation.
template class APrecond<AField<float,VECTOR> >;
template class APrecond<AField<double,VECTOR> >;

//============================================================END=====
